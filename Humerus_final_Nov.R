setwd("C:/Users/jcriosor/OneDrive/Dynamic Behaviour Lab/Hummingbirds/Geo_Morph_R/Hummer4")

### ============================================================
### BLOCK 1 — Base Data Import & Cleaning
### ============================================================

library(geomorph)
library(ape)
library(geiger)
library(dplyr)

### ------------------------------------------------------------
### 1. Load phylogeny
### ------------------------------------------------------------
tree <- read.tree("Trees24.tre")

# Standardize tree tip names to match metadata (spaces instead of underscores)
tree$tip.label <- gsub("_", " ", tree$tip.label)


### ------------------------------------------------------------
### 2. Load metadata
### ------------------------------------------------------------
metadata <- read.csv("metadata.csv", header = TRUE, sep = ";")
rownames(metadata) <- metadata$species


### ------------------------------------------------------------
### 3. Load TPS landmarks and assign species names
### ------------------------------------------------------------
shape_raw <- readland.tps(
  "humerus_male_RandL_cr3_forag_JCR.tps",
  specID = "imageID",
  warnmsg = FALSE
)

# Replace specimen IDs with species names from metadata
dimnames(shape_raw)[[3]] <- metadata$species


### ------------------------------------------------------------
### 4. Drop species in the tree that are not present in metadata
### ------------------------------------------------------------
tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(metadata)))


### ------------------------------------------------------------
### 5. Procrustes alignment + outlier removal
### ------------------------------------------------------------
gpa <- gpagen(
  shape_raw,
  PrinAxes = FALSE,
  ProcD = FALSE,
  print.progress = TRUE
)

shape_aligned <- gpa$coords

# Compute squared distances to the consensus (for outlier detection)
coords_mat <- two.d.array(shape_aligned)
mean_shape <- colMeans(coords_mat)

mean_matrix <- matrix(
  rep(mean_shape, each = nrow(coords_mat)),
  nrow = nrow(coords_mat),
  byrow = FALSE
)

distsq <- rowSums((coords_mat - mean_matrix)^2)

# Outlier filtering threshold (97.5%)
threshold <- quantile(distsq, 0.975)
keep_idx <- distsq < threshold

shape_clean  <- shape_aligned[, , keep_idx]
kept_ids     <- dimnames(shape_aligned)[[3]][keep_idx]
removed_ids  <- dimnames(shape_aligned)[[3]][!keep_idx]

cat("\nRemoved outliers:\n")
print(removed_ids)
cat("\nTotal specimens retained:", length(kept_ids), "\n")


### ------------------------------------------------------------
### 6. Prune tree and metadata to match shape_clean
### ------------------------------------------------------------
tree_clean <- keep.tip(tree, kept_ids)
metadata_clean <- metadata[kept_ids, ]

# Consistency checks
stopifnot(all(tree_clean$tip.label %in% rownames(metadata_clean)))
stopifnot(all(tree_clean$tip.label %in% dimnames(shape_clean)[[3]]))


### ------------------------------------------------------------
### 7. Compute centroid size (manual function)
### ------------------------------------------------------------
calc_centroid_size <- function(coords) {
  centroid <- colMeans(coords)
  diffs <- sweep(coords, 2, centroid)
  sqrt(sum(diffs^2))
}

centroid_sizes <- apply(shape_clean, 3, calc_centroid_size)
metadata_clean$centroid_size <- centroid_sizes


### ------------------------------------------------------------
### 8. Prepare ecological & morphological predictors
### ------------------------------------------------------------
bmass            <- metadata_clean$bmass
maxelev          <- as.numeric(metadata_clean$maxelev)
range_size       <- as.numeric(metadata_clean$Range.Size)
abs_centlatitude <- abs(as.numeric(metadata_clean$Centroid.Latitude))
foraging         <- as.factor(metadata_clean$foraging)

### ------------------------------------------------------------
### 9. Build geomorph data frame (gdf)
### ------------------------------------------------------------
gdf <- geomorph.data.frame(
  shape           = shape_clean,
  bmass           = bmass,
  maxelev         = maxelev,
  range           = range_size,
  abs_centlatitude = abs_centlatitude,
  centroid_size   = centroid_sizes,
  foraging        = foraging,
  tree            = tree_clean
)

### ============================================================
### END OF BLOCK 1
### ============================================================




### ============================================================
### BLOCK 2 – Multivariate PGLS of humeral shape (MAIN ANALYSIS)
### ============================================================

library(geomorph)
library(gtools)    # for permutations


### ------------------------------------------------------------
### 2.1 Candidate PGLS model set for humeral shape
### ------------------------------------------------------------

# All models use the same response (Procrustes shape) and phylogeny.
# Predictors are the ecological / morphological variables defined in gdf (Block 1).

model_set <- list(
  Minimal  = shape_clean ~ foraging,
  Flight   = shape_clean ~ foraging + maxelev + bmass,
  Extended = shape_clean ~ foraging + maxelev + bmass +
    centroid_size + range + abs_centlatitude
)

pgls_results <- data.frame(
  Model      = character(),
  Predictor  = character(),
  Rsq        = numeric(),
  p_value    = numeric(),
  stringsAsFactors = FALSE
)

for (model_name in names(model_set)) {
  form <- model_set[[model_name]]
  
  fit_try <- try(
    procD.pgls(
      f1    = form,
      data  = gdf,
      phy   = tree_clean,
      lambda = 0.99,          # fixed lambda used throughout the main analysis
      effect.type = "cohen",
      print.progress = FALSE
    ),
    silent = TRUE
  )
  
  if (inherits(fit_try, "try-error")) {
    cat("❌ Model failed:", model_name, "\n")
    next
  }
  
  aov_tab <- fit_try$aov.table
  
  # Drop "Residuals" and "Total" rows
  rows <- rownames(aov_tab)
  valid_rows <- which(!(rows %in% c("Residuals", "Total")))
  
  for (i in valid_rows) {
    pgls_results <- rbind(
      pgls_results,
      data.frame(
        Model     = model_name,
        Predictor = rownames(aov_tab)[i],
        Rsq       = round(aov_tab$Rsq[i], 4),
        p_value   = round(aov_tab$`Pr(>F)`[i], 4),
        stringsAsFactors = FALSE
      )
    )
  }
}

pgls_results


### ------------------------------------------------------------
### 2.2 Variable order optimisation via permutations
###     (sequential sums of squares / effect partitioning)
### ------------------------------------------------------------

# Here we explore how the order of predictors affects the sequential Rsq
# in procD.pgls. This is used as an exploratory tool; the final model is
# then fixed explicitly (see below) to keep the analysis reproducible.

predictors <- c(
  "foraging",
  "maxelev",
  "abs_centlatitude",
  "bmass",
  "centroid_size",
  "range"
)

perms <- permutations(
  n = length(predictors),
  r = length(predictors),
  v = predictors
)

perm_results_list <- list()

set.seed(999)  # reproducible permutations

for (i in 1:nrow(perms)) {
  
  this_order <- perms[i, ]
  formula_str <- paste(
    "shape_clean ~",
    paste(this_order, collapse = " + ")
  )
  
  fit <- procD.pgls(
    f1    = as.formula(formula_str),
    phy   = tree_clean,
    data  = gdf,
    iter  = 999,
    print.progress = FALSE
  )
  
  tab <- summary(fit)$table
  
  perm_results_list[[i]] <- list(
    order = this_order,
    R2    = tab$Rsq,
    p     = tab$`Pr(>F)`
  )
}

# Convert permutation results list into a data frame
perm_results_df <- do.call(
  rbind,
  lapply(perm_results_list, function(x) {
    data.frame(
      order = paste(x$order, collapse = " + "),
      R2_1  = x$R2[1],
      R2_2  = x$R2[2],
      R2_3  = x$R2[3],
      R2_4  = x$R2[4],
      R2_5  = x$R2[5],
      R2_6  = x$R2[6],
      p_1   = x$p[1],
      p_2   = x$p[2],
      p_3   = x$p[3],
      p_4   = x$p[4],
      p_5   = x$p[5],
      p_6   = x$p[6]
    )
  })
)

# Quick views (optional)
perm_results_df$R2_sum12 <- perm_results_df$R2_1 + perm_results_df$R2_2
perm_results_df$R2_var   <- apply(
  perm_results_df[, c("R2_1","R2_2","R2_3","R2_4","R2_5","R2_6")],
  1, var
)

best_by_p1   <- perm_results_df[order(perm_results_df$p_1), ]
best_by_R2   <- perm_results_df[order(-perm_results_df$R2_sum12), ]
best_stable  <- perm_results_df[order(perm_results_df$R2_var), ]

head(best_by_p1,   10)
head(best_by_R2,   10)
head(best_stable,  10)

# NOTE:
# Based on this exploration, we choose
# a biologically and statistically reasonable order for the final model:
# centroid_size + bmass + foraging + range + abs_centlatitude + maxelev


### ------------------------------------------------------------
### 2.3 Final multivariate PGLS model
### ------------------------------------------------------------

set.seed(2025)

final_formula <- shape_clean ~
  centroid_size +
  bmass +
  foraging +
  range +
  abs_centlatitude +
  maxelev

modelo_final <- procD.pgls(
  f1    = final_formula,
  data  = gdf,
  phy   = tree_clean,
  lambda = 0.99,         
  effect.type = "cohen",
  print.progress = FALSE
)

summary(modelo_final)



### ------------------------------------------------------------
### 2.5 Extract ANOVA table
### ------------------------------------------------------------

# Capture printed summary and isolate the ANOVA table lines
summary_capture <- capture.output(summary(modelo_final))

anova_lines <- summary_capture[
  grep("^\\s*(foraging|bmass|abs_centlatitude|centroid_size|range|maxelev)",
       summary_capture)
]

anova_table <- read.table(
  text  = anova_lines,
  header = FALSE,
  fill   = TRUE
)

# Assign column names manually (matching procD.pgls summary structure)
colnames(anova_table) <- c(
  "Predictor", "Df", "SS", "MS", "Rsq", "F", "Z", "p_value"
)

# Keep only the most relevant columns for the manuscript
anova_pub <- anova_table[, c("Predictor", "Df", "Rsq", "F", "Z", "p_value")]

print(anova_pub)

# Optional: export to CSV for direct use in the manuscript
write.csv(anova_pub,
          "statistics_final_PGLS_model_nov.csv",
          row.names = FALSE)

### ============================================================
### END OF BLOCK 2
### ============================================================



### ============================================================
### BLOCK 3 — RRphylo Evolutionary Rates of Humeral Shape
### ============================================================

library(RRphylo)
library(geomorph)
library(phytools)
library(dplyr)

### ------------------------------------------------------------
### 3.1 Prepare shape matrix ordered by tree tip labels
### ------------------------------------------------------------

# RRphylo requires a 2D matrix with rows = species, columns = coordinates
species_order <- tree_clean$tip.label

# Extract Procrustes shape array in the exact tip order
shape_array_ordered <- shape_clean[, , species_order]

# Convert to 2D matrix (landmarks × dims flattened)
Y_shape <- two.d.array(shape_array_ordered)
rownames(Y_shape) <- species_order

# Consistency check
stopifnot(all(rownames(Y_shape) == tree_clean$tip.label))


### ------------------------------------------------------------
### 3.2 Run RRphylo to estimate evolutionary rates
### ------------------------------------------------------------

RR_shape <- RRphylo(tree = tree_clean, y = Y_shape)

rates_all <- RR_shape$rates
ntip <- length(tree_clean$tip.label)

tip_branch_logical <- tree_clean$edge[, 2] <= ntip


MorphRate_RR <- RR_shape$rates[tip_branch_logical, 1]
names(MorphRate_RR) <- tree_clean$tip.label

cat("Total species with rate estimates:", length(MorphRate_RR), "\n")
summary(MorphRate_RR)


### ------------------------------------------------------------
### 3.3 Visualisation of rate distribution
### ------------------------------------------------------------

library(ggplot2)

df_rates <- data.frame(
  species = names(MorphRate_RR),
  rate    = MorphRate_RR,
  lograte = log10(MorphRate_RR)
)

ggplot(df_rates, aes(x = lograte)) +
  geom_histogram(aes(y = ..density..),
                 bins = 25, color = "black", fill = "grey80", alpha = 0.8) +
  geom_density(size = 1.2) +
  theme_classic(base_size = 14) +
  labs(
    title = "Distribution of evolutionary rates (RRphylo)",
    x = "log10 morphological rate",
    y = "Density"
  )


### ------------------------------------------------------------
### 3.4 filter extreme outliers for regression purposes
### ------------------------------------------------------------

# Conservative filter: remove 1.25% highest and lowest extremes
q_low  <- quantile(MorphRate_RR, 0.0125, na.rm = TRUE)
q_high <- quantile(MorphRate_RR, 0.9875, na.rm = TRUE)

valid_rates <- MorphRate_RR >= q_low & MorphRate_RR <= q_high
MorphRate_RR_filt <- MorphRate_RR[valid_rates]

cat("Rates retained after filtering:", sum(valid_rates), "of", length(MorphRate_RR), "\n")


### ------------------------------------------------------------
### 3.5 Identify outlier species (IQR-based)
### ------------------------------------------------------------

Q1 <- quantile(MorphRate_RR, 0.25)
Q3 <- quantile(MorphRate_RR, 0.75)
IQR_val <- Q3 - Q1

upper_lim <- Q3 + 1.5 * IQR_val
lower_lim <- Q1 - 1.5 * IQR_val

outliers_high <- MorphRate_RR[MorphRate_RR > upper_lim]
outliers_low  <- MorphRate_RR[MorphRate_RR < lower_lim]

cat("High-rate outliers:\n")
print(outliers_high)
cat("Low-rate outliers:\n")
print(outliers_low)


### ------------------------------------------------------------
### 3.6 Build species-level dataframe for rate analyses
### ------------------------------------------------------------

species_df <- data.frame(
  species         = tree_clean$tip.label,
  MorphRate       = MorphRate_RR,
  bmass           = gdf$bmass,
  range           = gdf$range,
  centroid_size   = gdf$centroid_size,
  maxelev         = gdf$maxelev,
  abs_centlatitude = gdf$abs_centlatitude,
  foraging        = gdf$foraging,
  stringsAsFactors = FALSE
)

# Filtered version (for regression models)
species_df_filt <- species_df[valid_rates, ]

# Quick inspection
str(species_df_filt)
summary(species_df_filt$MorphRate)


### ------------------------------------------------------------
### 3.7 Phylogenetic signal in evolutionary rates
### ------------------------------------------------------------

lambda_rate <- phylosig(
  tree_clean,
  MorphRate_RR,
  method = "lambda",
  test = TRUE
)

cat("\n=== Phylogenetic signal of evolutionary rates ===\n")
print(lambda_rate)

### ============================================================
### END OF BLOCK 3
### ============================================================



### ============================================================
### BLOCK 4 – Evolutionary rate correlates + group-wise rates
### ============================================================

library(phylolm)
library(geomorph)
library(reshape2)
library(dplyr)

### ------------------------------------------------------------
### 4.1 – Univariate PGLS models: MorphRate ~ predictors
###       (using RRphylo rates)
### ------------------------------------------------------------

# Continuous predictors to be tested against RRphylo rates
rate_predictors <- c("bmass", "centroid_size", "range",
                     "maxelev", "abs_centlatitude")

rate_models_list <- list()
rate_pgls_table  <- data.frame()

for (pred in rate_predictors) {
  
  # Build formula MorphRate ~ pred
  form <- as.formula(paste("MorphRate ~", pred))
  
  # Subset filtered dataframe to the variables needed
  df_mod <- species_df_filt[, c("MorphRate", pred, "species")]
  
  # Match phylogeny to the species present in df_mod
  phy_sub <- drop.tip(
    tree_clean,
    setdiff(tree_clean$tip.label, df_mod$species)
  )
  
  # Fit PGLS with lambda estimated by ML
  m <- phylolm(
    form,
    data = df_mod,
    phy  = phy_sub,
    model = "lambda"
  )
  
  rate_models_list[[pred]] <- m
  
  # Extract statistics for the predictor
  s <- summary(m)
  coef_row <- s$coefficients[2, ]  # second row = focal predictor
  
  coef_est <- coef_row["Estimate"]
  coef_se  <- coef_row["StdErr"]
  coef_p   <- coef_row["p.value"]
  lambda_m <- m$optpar  # ML lambda
  
  rate_pgls_table <- rbind(
    rate_pgls_table,
    data.frame(
      Predictor = pred,
      Estimate  = coef_est,
      StdError  = coef_se,
      p_value   = coef_p,
      Lambda    = lambda_m,
      stringsAsFactors = FALSE
    )
  )
}

rate_pgls_table


### ------------------------------------------------------------
### 4.2 – Final table of PGLS (rates vs predictors) 
### ------------------------------------------------------------

table_rates_final <- rate_pgls_table %>%
  rename(
    `Std. Error` = StdError,
    `p-value`    = p_value,
    `Lambda (ML)` = Lambda
  )

print(table_rates_final)

# Optional export
write.csv(
  table_rates_final,
  "Table_Rates_PGLS_vs_Traits.csv",
  row.names = FALSE
)


### ------------------------------------------------------------
### 4.3 – Evolutionary rates by foraging strategy (compare.evol.rates)
###       MAIN GROUP-WISE RATE ANALYSIS
### ------------------------------------------------------------

# Here we test whether humeral shape evolves at different rates among
# foraging regimes (ter, tra, opp). We use the same Procrustes shape
# array and phylogeny as in previous blocks.

# Ensure group factor has names matching the specimens in shape_clean
foraging_named <- foraging
names(foraging_named) <- dimnames(shape_clean)[[3]]


ER_foraging <- compare.evol.rates(
  A      = shape_clean,
  phy    = tree_clean,
  gp     = foraging_named,
  iter   = 999,
  method = "permutation",
  print.progress = FALSE
)


ER_foraging   # full output


# Extract group-specific sigmas
sigmas <- ER_foraging$sigma.d.gp

# Pairwise rate ratios and p-values
# (ER_foraging$rate.ratio and ER_foraging$P.value are usually given
# as vectors; we reorganise them into 3x3 matrices for clarity.)

groups <- names(sigmas)

# Initialise matrices
ratios_mat <- matrix(NA, nrow = length(groups), ncol = length(groups),
                     dimnames = list(groups, groups))

pvals_mat <- matrix(
  ER_foraging$P.value,
  nrow = length(groups), ncol = length(groups),
  dimnames = list(groups, groups)
)

# Fill the ratio matrix as sigma_i / sigma_j
for (i in seq_along(groups)) {
  for (j in seq_along(groups)) {
    ratios_mat[i, j] <- sigmas[groups[i]] / sigmas[groups[j]]
  }
}

# Convert both matrices to long format
ratios_long <- melt(
  ratios_mat,
  varnames = c("Group1", "Group2"),
  value.name = "Rate_Ratio"
)

pvals_long <- melt(
  pvals_mat,
  varnames = c("Group1", "Group2"),
  value.name = "p_value"
)

# Merge into a single table
table_ER_foraging <- merge(
  ratios_long,
  pvals_long,
  by = c("Group1", "Group2")
)

# Remove diagonal (comparisons of a group with itself)
table_ER_foraging <- table_ER_foraging[
  table_ER_foraging$Group1 != table_ER_foraging$Group2, ]

# Order for readability
table_ER_foraging <- table_ER_foraging[
  order(table_ER_foraging$Group1, table_ER_foraging$Group2), ]

table_ER_foraging

# Optional export
write.csv(
  table_ER_foraging,
  "Table_ER_Foraging_RateRatios.csv",
  row.names = FALSE
)

### ============================================================
### END OF BLOCK 4
### ============================================================




### ============================================================
### BLOCK 5 — Manuscript Figures (Shape Space & Integration)
### ============================================================

library(geomorph)
library(ggplot2)
library(phytools)
library(dplyr)

paleta_foraging <- c(
  "ter" = "#440154",   # morado viridis
  "tra" = "#21908C",   # verde viridis
  "opp" = "#F46D43"    # naranja intenso
)

### ------------------------------------------------------------
### 5.1 — PCA of the final multivariate PGLS model (shape space)
### ------------------------------------------------------------

# This generates a PCA of fitted values from the final PGLS model.
pc_plot <- plot(
  modelo_final,
  type = "PC",
  pch = 21,
  bg = paleta_foraging[gdf$foraging],
  cex = 1.2
)

# Add convex hulls (one per foraging group)
shapeHulls(
  pc_plot,
  groups = gdf$foraging,
  group.cols = paleta_foraging,
  group.lwd = c(2,2,2),
  group.lty = c(1,1,1)
)

legend(
  "topright",
  legend = levels(gdf$foraging),
  col = paleta_foraging,
  pch = 19,
  pt.cex = 1.3,
  bty = "n"
)


### ------------------------------------------------------------
### 5.2 — Phylomorphospace (GM + phylogeny)
### ------------------------------------------------------------

gm_res <- gm.prcomp(shape_clean, phy = tree_clean)

# Align foraging vector to specimen order
specimen_names <- dimnames(shape_clean)[[3]]
foraging_plot <- metadata_clean[specimen_names, "foraging"]

# Convert to factor with correct levels
foraging_plot <- factor(foraging_plot, levels = levels(gdf$foraging))


plot(
  gm_res,
  phylo = TRUE,
  pch = 21,
  cex = 1.5,
  bg = paleta_foraging[foraging_plot],
  phylo.par = list(
    tip.labels = FALSE,
    node.labels = FALSE,
    node.cex = 0.2
  )
)

legend(
  "topright",
  legend = levels(foraging_plot),
  col = paleta_foraging,
  pch = 19,
  pt.cex = 1.3,
  bty = "n"
)



### ------------------------------------------------------------
### 5.3 — Phylogenetic signal in humeral shape (K + plot)
### ------------------------------------------------------------

PS_shape <- physignal(
  gdf$shape,
  gdf$tree,
  print.progress = FALSE
)

summary(PS_shape)
plot(PS_shape)


### ------------------------------------------------------------
### 5.4 — Phenotypic integration: Proximal vs. Distal regions (PLS)
### ------------------------------------------------------------

# Partition landmarks into two blocks:
mean_shape <- mshape(shape_clean)
plot(mean_shape, pch=16)
text(mean_shape, labels = 1:nrow(mean_shape), pos=3)


# Replace the sequence below with the TRUE proximal/distal assignment
proximal_LMs <- c(1,2,3,4,5,6,7,8,9)      # EJEMPLO, AJUSTAR
distal_LMs   <- c(10,11,12,13,14,15,16,
                  17,18,19,20,21,22,23,24,25)   # EJEMPLO, AJUSTAR

land.gps <- rep(NA, nrow(mean_shape))
land.gps[proximal_LMs] <- "A"
land.gps[distal_LMs]   <- "B"


if (any(is.na(land.gps)))
  stop("Some landmarks remain unassigned. Fix proximal/distal lists.")

pls_res <- phylo.integration(
  A = gdf$shape,
  partition.gp = land.gps,
  phy = tree_clean,
  print.progress = FALSE
)

summary(pls_res)

scores_block1 <- pls_res$XScores  # proximal scores
scores_block2 <- pls_res$YScores  # distal scores

# Build dataframe
pls_df <- data.frame(
  Block1 = scores_block1[, 1],
  Block2 = scores_block2[, 1]
)

# Remove extreme points (> ±3 SD)
pls_df_clean <- pls_df %>%
  mutate(
    z1 = scale(Block1),
    z2 = scale(Block2)
  ) %>%
  filter(abs(z1) < 3, abs(z2) < 3)

# Plot integration
ggplot(pls_df_clean, aes(x = Block1, y = Block2)) +
  geom_point(size = 2.5, alpha = 0.7, color = "steelblue") +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "purple",
    size = 1.2
  ) +
  theme_minimal(base_size = 16) +
  labs(
    title = "Phenotypic integration between proximal and distal humeral regions",
    x = "PLS1 — Proximal region",
    y = "PLS1 — Distal region"
  )


### ============================================================
### END OF BLOCK 5
### ============================================================




### ============================================================
### 5.1 – Allometric correction for CVA
### ============================================================

library(geomorph)

# Allometric object
log_centroid <- log(centroid_sizes)

# Allometric model
fit_allom <- procD.lm(shape_clean ~ log_centroid, iter = 0)

# Allometric residuals
shape_allom_res <- arrayspecs(residuals(fit_allom), 
                              p = dim(shape_clean)[1], 
                              k = dim(shape_clean)[2])




### ============================================================
### 5.2 – CVA using leave-one-out cross-validation
### ============================================================

library(geomorph)
library(Morpho)

# 1. Allometry
log_centroid <- log(centroid_sizes)
fit_allom <- procD.lm(shape_clean ~ log_centroid, iter = 0)

shape_allom_res <- arrayspecs(
  residuals(fit_allom),
  p = dim(shape_clean)[1],
  k = dim(shape_clean)[2]
)

# 2. factor names
names(foraging) <- dimnames(shape_clean)[[3]]

# 3. CVA (Morpho)
library(Morpho)
CVA_out <- CVA(shape_allom_res,
               group = foraging,
               plot = FALSE,
               cv = TRUE)

CVA_out

### ============================================================
### 5.3 – Effect of Foraging regime in shape
### ============================================================

foraging_model <- procD.lm(shape_allom_res ~ foraging, 
                           iter = 999, 
                           print.progress = FALSE)

summary(foraging_model)

### ============================================================
### END OF BLOCK 5.2
### ============================================================




### ============================================================
### BLOCK 6 — Partial Effect Plots (PGLS model effects)
### ============================================================

library(ggplot2)
library(dplyr)
library(patchwork)

### ------------------------------------------------------------
### 6.1 — Project fitted values into 1D morphospace (PC1)
### ------------------------------------------------------------

# Fitted values from the final PGLS model
fitted_vals <- modelo_final$pgls.fitted

# PCA on fitted multivariate predictions
pca_fit <- prcomp(fitted_vals, scale. = TRUE)
PC1 <- pca_fit$x[, 1]

# Build dataframe with predictors + PC1
plot_df <- as.data.frame(modelo_final$data)
plot_df$Morpho_PC1 <- PC1


### ------------------------------------------------------------
### 6.2 — Define predictors to plot (same as model_final)
### ------------------------------------------------------------

vars_to_plot <- c(
  "foraging",
  "bmass",
  "centroid_size",
  "range",
  "maxelev",
  "abs_centlatitude"
)

# Colors for categorical variables (from global palette)
paleta_extendida <- paleta_foraging


### ------------------------------------------------------------
### 6.3 — Function to generate each panel
### ------------------------------------------------------------

partial_plot <- function(df, var) {
  
  # Remove extreme outliers for cleaner visualisation
  df_clean <- df
  
  if (var == "range") {
    df_clean <- df_clean %>% filter(range <= quantile(range, 0.95, na.rm = TRUE))
  }
  
  if (var == "centroid_size") {
    df_clean <- df_clean %>% filter(centroid_size <= quantile(centroid_size, 0.99, na.rm = TRUE))
  }
  
  if (var == "abs_centlatitude") {
    df_clean <- df_clean %>% filter(abs_centlatitude <= quantile(abs_centlatitude, 0.99, na.rm = TRUE))
  }
  
  # Categorical predictors (foraging)
  if (is.factor(df_clean[[var]]) || is.character(df_clean[[var]])) {
    
    return(
      ggplot(df_clean, aes_string(x = var, y = "Morpho_PC1", fill = var, color = var)) +
        geom_boxplot(alpha = 0.6, outlier.shape = NA) +
        geom_jitter(position = position_jitter(width = 0.2), size = 2.2, alpha = 0.6) +
        scale_fill_manual(values = paleta_extendida) +
        scale_color_manual(values = paleta_extendida) +
        theme_minimal(base_size = 15) +
        labs(
          title = paste("Morphospace vs", var),
          x = var,
          y = "Humeral shape (PC1)"
        )
    )
    
  } else {
    # Continuous predictors
    return(
      ggplot(df_clean, aes_string(x = var, y = "Morpho_PC1")) +
        geom_point(size = 2.2, alpha = 0.6, color = "steelblue") +
        geom_smooth(method = "lm", se = TRUE, color = "purple", size = 1.2) +
        theme_minimal(base_size = 15) +
        labs(
          title = paste("Morphospace vs", var),
          x = var,
          y = "Humeral shape (PC1)"
        )
    )
  }
}


### ------------------------------------------------------------
### 6.4 — Generate all panels and assemble figure
### ------------------------------------------------------------

plot_list <- lapply(vars_to_plot, function(v) partial_plot(plot_df, v))
names(plot_list) <- vars_to_plot

# Combine panels using patchwork
full_partial_plot <- wrap_plots(plot_list, ncol = 2)

full_partial_plot

### ============================================================
### END OF BLOCK 6
### ============================================================



### ============================================================
### BLOCK 6.2 — Figures for RRphylo Rates
### ============================================================

library(ggplot2)
library(dplyr)
library(patchwork)

### ------------------------------------------------------------
### Build unified dataframe for rate plots
### ------------------------------------------------------------

df_rateplot <- species_df %>%
  mutate(
    lograte = log10(MorphRate),
    foraging = factor(foraging, levels = levels(gdf$foraging))
  )


### ------------------------------------------------------------
### 6.2.1 — Distribution of RRphylo rates
### ------------------------------------------------------------

p_rate_dist <- ggplot(df_rateplot, aes(x = lograte)) +
  geom_histogram(aes(y = ..density..),
                 bins = 25,
                 fill = "grey85",
                 color = "black") +
  geom_density(size = 1.2, color = "purple") +
  theme_classic(base_size = 16) +
  labs(
    title = "Distribution of evolutionary rates (RRphylo)",
    x = "log10 morphological rate",
    y = "Density"
  )


### ------------------------------------------------------------
### 6.2.2 — Rates by foraging strategy
### ------------------------------------------------------------

p_rate_foraging <- ggplot(df_rateplot,
                          aes(x = foraging, y = lograte,
                              fill = foraging, color = foraging)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  scale_fill_manual(values = paleta_foraging) +
  scale_color_manual(values = paleta_foraging) +
  theme_minimal(base_size = 16) +
  labs(
    title = "Evolutionary rates by foraging strategy",
    x = "Foraging regime",
    y = "log10 rate"
  )


### ------------------------------------------------------------
### 6.2.3 — Partial effect plots (Rate ~ traits)
### ------------------------------------------------------------

rate_vars <- c("bmass", "centroid_size", "range",
               "maxelev", "abs_centlatitude")

partial_rate_plot <- function(df, var) {
  
  # For aesthetics: remove extreme values only in the predictor
  df_clean <- df %>%
    filter((!!as.name(var)) <= quantile((!!as.name(var)), 0.99, na.rm = TRUE))
  
  ggplot(df_clean, aes_string(x = var, y = "lograte")) +
    geom_point(size = 2.2, alpha = 0.6, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "purple", size = 1.2) +
    theme_minimal(base_size = 15) +
    labs(
      title = paste("Evolutionary rate vs", var),
      x = var,
      y = "log10 RRphylo rate"
    )
}

rate_panels <- lapply(rate_vars, function(v) partial_rate_plot(df_rateplot, v))
names(rate_panels) <- rate_vars

p_rate_effects <- wrap_plots(rate_panels, ncol = 2)


### ------------------------------------------------------------
### 6.2.4 — Phylogenetic signal in RRphylo rates
### ------------------------------------------------------------

lambda_rate <- phylosig(tree_clean, MorphRate_RR, method = "lambda", test = TRUE)

sim_null <- replicate(
  999,
  phylosig(tree_clean, sample(MorphRate_RR), method = "lambda", test = FALSE)
)

### Safe simulation of lambda values
sim_vals <- replicate(999, {
  out <- try(phylosig(tree_clean,
                      sample(MorphRate_RR),
                      method = "lambda",
                      test = FALSE),
             silent = TRUE)
  
  # if phylosig failed → return NA
  if (inherits(out, "try-error")) return(NA)
  
  # if phylosig returned a list (common)
  if (is.list(out)) {
    if ("lambda" %in% names(out)) return(out$lambda)
    return(NA)
  }
  
  # if phylosig returned a numeric value
  if (is.numeric(out)) return(out)
  
  # fallback
  return(NA)
})

df_signal <- data.frame(
  lambda = c(sim_vals, lambda_rate$lambda),
  type   = c(rep("null", length(sim_vals)), "observed")
)



p_rate_signal <- ggplot(df_signal, aes(x = lambda)) +
  geom_histogram(data = subset(df_signal, type == "null"),
                 bins = 30,
                 fill = "grey85",
                 color = "black") +
  geom_vline(xintercept = lambda_rate$lambda,
             color = "red",
             size = 1.2) +
  theme_classic(base_size = 16) +
  labs(
    title = "Phylogenetic signal in evolutionary rates",
    subtitle = paste("Observed λ =", round(lambda_rate$lambda, 3)),
    x = "Lambda",
    y = "Frequency"
  )


### ------------------------------------------------------------
### Show all figures
### ------------------------------------------------------------

p_rate_dist
p_rate_foraging
p_rate_effects
p_rate_signal

### ============================================================
### END OF BLOCK 6.2
### ============================================================


