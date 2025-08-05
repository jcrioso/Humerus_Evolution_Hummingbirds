library(geomorph)
library(geiger)
library(ape)
library(ggplot2)
library(phytools)

# Paleta definitiva para grupos de forrajeo
paleta_foraging <- c(
  "ter" = "#440154",   # morado viridis
  "tra" = "#21908C",   # verde viridis
  "opp" = "#F46D43"    # naranja intenso
)


# Cargar árbol
tree <- read.tree("Trees24.tre")
plot(tree, 
     type = "fan", 
     cex = 0.6, 
     label.offset = 0.01,
     no.margin = TRUE,
     edge.color = "darkgrey",
     tip.color = "black",
     edge.width = 1.2)


# Convertir nombres del árbol al mismo formato que en pcs$species
tree$tip.label <- gsub("_", " ", tree$tip.label)

# Cargar metadata
metadata <- read.csv("metadata.csv", header = TRUE, sep = ";")

# Cargar landmarks desde archivo TPS
shape <- readland.tps("humerus_male_RandL_cr3_forag_JCR.tps", specID = "imageID", warnmsg = FALSE)


# Renombrar los especímenes del array shape con species
dimnames(shape)[[3]] <- metadata$species

rownames(metadata) <- metadata$species


# Crear geomorph data frame alineado con el árbol
tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(metadata)))
# Definir conectividad entre landmarks según el outline
links <- matrix(c(
  1, 2,
  2, 3,
  3, 4,
  4, 5,
  5, 6,
  6, 7,
  7, 8,
  8, 9,
  9, 10,
  10, 11,
  11, 12,
  12, 13,
  13, 14,
  14, 15,
  15, 16,
  16, 17,
  17, 18,
  18, 19,
  19, 20,
  20, 21,
  21, 22,
  22, 23,
  23, 24,
  24, 25,
  25, 1,    
  13, 16
), ncol = 2, byrow = TRUE)


#### procrustes
gpa <- gpagen(shape, PrinAxes = FALSE, ProcD = FALSE, print.progress = TRUE)
shape_aligned <- gpa$coords

# Calcular distancias cuadradas al consenso
coords_mat <- two.d.array(shape_aligned)
mean_shape <- colMeans(coords_mat)
# Asegurar que mean_shape tenga misma forma que coords_mat
mean_matrix <- matrix(rep(mean_shape, each = nrow(coords_mat)), nrow = nrow(coords_mat), byrow = FALSE)

# Calcular distancia euclidiana cuadrada por fila (especímen)
distsq <- rowSums((coords_mat - mean_matrix)^2)

# Eliminar outliers: conservar el 97.5% más cercanos al consenso
threshold <- quantile(distsq, 0.975)
keep <- distsq < threshold


specimen_ids <- dimnames(shape_aligned)[[3]]
kept_ids <- specimen_ids[keep]
removed_ids <- specimen_ids[!keep]

# Crear shape limpio
shape_clean <- shape_aligned[, , keep]

cat("✅ Se eliminaron", length(removed_ids), "especímenes como outliers (2.5% más extremos)\n\n")
cat("🔴 Especímenes eliminados:\n")
print(removed_ids)

cat("\n✅ Se conservaron", length(kept_ids), "especímenes\n")

dev.off()

# Visualización limpia
plotAllSpecimens(shape_clean, links = links, label = TRUE)
plotAllSpecimens(shape_clean, label = FALSE)



#### Outline adicional

plot(NA, xlim = range(shape_clean[, 1, ]), ylim = range(shape_clean[, 2, ]),
     asp = 1, xlab = "", ylab = "", axes = FALSE, main = "Forma promedio y dispersión de landmarks")

for (i in 1:dim(shape_clean)[3]) {
  specimen <- shape_clean[, , i]
  points(specimen, pch = 16, col = rgb(0.3, 0.3, 0.3, 0.9), cex = 0.2)  # azul marino translúcido
}

mean_shape <- mshape(shape_clean)

for (j in 1:nrow(links)) {
  segments(mean_shape[links[j, 1], 1], mean_shape[links[j, 1], 2],
           mean_shape[links[j, 2], 1], mean_shape[links[j, 2], 2],
           col = "navy", lwd = 2)
}
points(mean_shape, col = "red", pch = 16, cex = 0.5)




##### Outline adicional 2

plot(NA, xlim = range(shape_clean[, 1, ]), ylim = range(shape_clean[, 2, ]),
     asp = 1, xlab = "", ylab = "", axes = FALSE, main = "Outline de todos los especímenes")

for (i in 1:dim(shape_clean)[3]) {
  specimen <- shape_clean[, , i]
  
  for (j in 1:nrow(links)) {
    segments(specimen[links[j, 1], 1], specimen[links[j, 1], 2],
             specimen[links[j, 2], 1], specimen[links[j, 2], 2],
             col = rgb(0, 0, 0.5, 0.1), lwd = 1)
  }
  
  
  points(specimen, pch = 16, col = rgb(0, 0, 0.5, 0.1), cex = 0.3)
}

mean_shape <- mshape(shape_clean)

for (j in 1:nrow(links)) {
  segments(mean_shape[links[j, 1], 1], mean_shape[links[j, 1], 2],
           mean_shape[links[j, 2], 1], mean_shape[links[j, 2], 2],
           col = "navy", lwd = 2)
}

points(mean_shape, col = "red", pch = 16, cex = 0.5)


## alinear arbol con metadatos sin outliers
tree_clean <- keep.tip(tree, kept_ids)

cat("🧩 Árbol ajustado tiene", length(tree_clean$tip.label), "especies.\n")
cat("✅ Shape limpio tiene", dim(shape_clean)[3], "especies.\n")

plot(tree_clean,
     main = "Árbol con especies conservadas",
     type = "fan", 
     cex = 0.6, 
     label.offset = 0.01,
     no.margin = TRUE,
     edge.color = "darkgrey",
     tip.color = "black",
     edge.width = 1.2)

metadata_clean <- metadata[kept_ids, ]



##### calcular elevación como variable binaria y trinomial para algunos analisis

low_perc <- quantile(metadata_clean$maxelev, 0.33, na.rm = TRUE)
high_perc <- quantile(metadata_clean$maxelev, 0.66, na.rm = TRUE)

metadata_clean$elev_bin <- ifelse(metadata_clean$maxelev <= high_perc, "low", "high")
metadata_clean$elev_bin <- as.factor(metadata_clean$elev_bin)

metadata_clean$elev_trinomial <- cut(
  metadata_clean$maxelev,
  breaks = c(-Inf, low_perc, high_perc, Inf),
  labels = c("low", "medium", "high")
)
metadata_clean$elev_trinomial <- as.factor(metadata_clean$elev_trinomial)

elev_bin <- metadata_clean$elev_bin
names(elev_bin) <- rownames(metadata_clean)
elev_bin <- elev_bin[names(elev_bin) %in% dimnames(shape_clean)[[3]]]
elev_bin <- elev_bin[match(dimnames(shape_clean)[[3]], names(elev_bin))]

elev_tri <- metadata_clean$elev_trinomial
names(elev_tri) <- rownames(metadata_clean)
elev_tri <- elev_tri[names(elev_tri) %in% dimnames(shape_clean)[[3]]]
elev_tri <- elev_tri[match(dimnames(shape_clean)[[3]], names(elev_tri))]

## convertir a latitud a absoluta para usarla en todos los modelos y figuras
metadata_clean$abs_centlatitude <- abs(metadata_clean$Centroid.Latitude)

# Extraer variables limpias
bmass <- metadata_clean$bmass
avelev <- metadata_clean$avelev
maxelev <- as.integer(metadata_clean$maxelev)
foraging <- as.factor(metadata_clean$foraging)
clade <- as.factor(metadata_clean$clade)
maxlatitude <- as.numeric(metadata_clean$Max.Latitude)
range <- as.numeric(metadata_clean$Range.Size)
keel <- as.numeric(metadata_clean$keel)
centlatitude <- as.numeric(metadata_clean$Centroid.Latitude)
abs_centlatitude <- as.numeric(metadata_clean$abs_centlatitude)
longitude <- as.numeric(metadata_clean$Centroid.Longitude)
habitat <- as.factor(metadata_clean$habitat)
habitat_density <- as.factor(metadata_clean$habitat.density)
elev_bin <- as.factor(metadata_clean$elev_bin)
elev_tri <- as.factor(metadata_clean$elev_trinomial)

gdf <- geomorph.data.frame(
  shape = shape_clean,
  bmass = bmass,
  avelev = avelev,
  maxelev = maxelev,
  elev_bin = elev_bin,
  elev_tri = elev_tri,
  foraging = foraging,
  clade = clade,
  maxlatitude = maxlatitude,
  range = range,
  keel = keel,
  centlatitude = centlatitude,
  abs_centlatitude = abs_centlatitude,
  longitude = longitude,
  habitat = habitat,
  habitat_density = habitat_density,
  tree = tree_clean
)


# Función para calcular manualmente el centroid size
calc_centroid_size <- function(coords) {
  centroid <- colMeans(coords)
  diffs <- sweep(coords, 2, centroid)
  sqrt(sum(diffs^2))
}

centroid_sizes <- apply(gdf$shape, 3, calc_centroid_size)

metadata_clean$centroid_size <- centroid_sizes
head(metadata_clean$centroid_size)


centroid_size <- as.numeric((metadata_clean$centroid_size))

gdf <- geomorph.data.frame(
  shape = shape_clean,
  bmass = bmass,
  avelev = avelev,
  maxelev = maxelev,
  elev_bin = elev_bin,
  elev_tri = elev_tri,
  foraging = foraging,
  clade = clade,
  maxlatitude = maxlatitude,
  range = range,
  keel = keel,
  centlatitude = centlatitude,
  abs_centlatitude = abs_centlatitude,
  longitude = longitude,
  centroid_size = centroid_size,
  habitat = habitat,
  habitat_density = habitat_density,
  tree = tree_clean
)



# PGLS

# matriz de correlaciones
num_vars_updated <- metadata_clean[, c("maxelev", "bmass", "Range.Size", 
                                       "Centroid.Longitude", "abs_centlatitude", 
                                       "centroid_size")]

cor_matrix_updated <- cor(num_vars_updated, use = "complete.obs")
print(round(cor_matrix_updated, 3))
library(corrplot)
corrplot(cor_matrix_updated, method = "color", type = "upper", addCoef.col = "black", tl.cex = 0.8)

## el centroide de la latitud esta muy correlacionado con el centride de la longitud, 
## indicios de quedarnos con unos solo. Vamos a quedarnos con Latitud


# VIF
library(car)
library(ggplot2)
library(dplyr)
modelo_vif_updated <- lm(centroid_size ~ foraging + maxelev + 
                           + abs_centlatitude + Range.Size + habitat + 
                           habitat.density, data = metadata_clean)
vif_valores_updated <- vif(modelo_vif_updated)
print(vif_valores_updated)


#### Selección de modelos candidatos
modelos <- list(
  "Minimo" = shape_clean ~ foraging,
  "Vuelo" = shape_clean ~ foraging + maxelev + bmass,
  "Extendido" = shape_clean ~ foraging + maxelev + bmass + centroid_size + range + abs_centlatitude,
  "Completo" = shape_clean ~ foraging + maxelev + bmass + centroid_size + range + abs_centlatitude + habitat + habitat_density
)

resultados <- data.frame(
  Modelo = character(),
  Predictora = character(),
  Rsq = numeric(),
  P_valor = numeric(),
  stringsAsFactors = FALSE
)

# Loop por cada modelo en la lista 'modelos'
for (nombre_modelo in names(modelos)) {
  formula_actual <- modelos[[nombre_modelo]]
  
  modelo <- try(procD.pgls(f1 = formula_actual,
                           effect.type = 'cohen',
                           data = gdf,
                           phy = tree_clean,
                           lambda = 0.98,
                           print.progress = FALSE),
                silent = TRUE)
  
  if (inherits(modelo, "try-error")) {
    cat("❌ Modelo fallido:", nombre_modelo, "\n")
    next
  }
  
  aov_table <- modelo$aov.table
  
  for (i in 1:(nrow(aov_table) - 2)) {
    predictor_name <- rownames(aov_table)[i]
    rsq_val <- aov_table$Rsq[i]
    p_val <- aov_table$`Pr(>F)`[i]
    
    resultados <- rbind(resultados, data.frame(
      Modelo = nombre_modelo,
      Predictora = predictor_name,
      Rsq = round(rsq_val, 4),
      P_valor = round(p_val, 4),
      stringsAsFactors = FALSE
    ))
  }
}

print(resultados)




####### mejor modelo

modelo_final <- procD.pgls(f1 = shape_clean ~ foraging + abs_centlatitude + bmass + centroid_size + range + maxelev,
                           effect.type = 'cohen',
                           data = gdf,
                           phy = tree_clean,
                           lambda = 0.99,
                           print.progress = FALSE)
summary(modelo_final)


##### verificación de lambda ideal

logLik_lambda <- function(lambda_value) {
  modelo <- try(procD.pgls(
    f1 = shape_clean ~ foraging + maxelev + abs_centlatitude + bmass + centroid_size + Range.Size,
    effect.type = 'cohen',
    data = gdf,
    phy = tree_clean,
    lambda = lambda_value,
    print.progress = FALSE
  ), silent = TRUE)
  
  if (inherits(modelo, "try-error")) {
    return(Inf)
  }
  
  # procD.pgls no da logLik directamente, pero podemos usar AIC para derivarlo:
  # AIC = -2 * logLik + 2 * k
  # Despejamos logLik:
  aic_val <- modelo$aic
  
  k <- length(coef(modelo))
  
  logLik_val <- -(aic_val - 2 * k) / 2
  
  return(-logLik_val)
}

lambda_opt <- optimize(logLik_lambda, interval = c(0, 1))

cat("✅ Lambda óptimo:", lambda_opt$minimum, "\n")
cat("🔍 Log-verosimilitud (negativa):", lambda_opt$objective, "\n") ### Lambda optimo = 0.99



# tabular

captura <- capture.output(summary(modelo_final))
lineas_tabla <- captura[grep("^\\s*(foraging|bmass|abs_centlatitude|centroid_size|range|habitat|maxelev)", captura)]

tabla_publicacion <- read.table(text = lineas_tabla, header = FALSE, fill = TRUE)

colnames(tabla_publicacion) <- c("Predictora", "Df", "SS", "MS", "Rsq", "F", "Z", "p-valor")

tabla_publicacion <- tabla_publicacion[, c("Predictora", "Df", "Rsq", "F", "Z", "p-valor")]
tabla_publicacion <- tabla_publicacion[-c(7:8),]
print(tabla_publicacion)
write.csv(tabla_publicacion, "statistics_modelo_final.csv")




### graficar tamaño del efecto y el valor de p

tabla_publicacion$Rsq <- as.numeric(as.character(tabla_publicacion$Rsq))
tabla_publicacion$`p-valor` <- as.numeric(as.character(tabla_publicacion$`p-valor`))

ggplot(tabla_publicacion, aes(x = reorder(Predictora, Rsq), y = Rsq)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Tamaño del efecto (R²) por predictor",
    x = "Predictor",
    y = "R² explicado"
  ) +
  geom_text(aes(label = round(Rsq, 3)), hjust = -0.1, size = 4) +
  ylim(0, max(tabla_publicacion$Rsq, na.rm = TRUE) * 1.2)

ggplot(tabla_publicacion, aes(x = reorder(Predictora, -`p-valor`), y = `p-valor`)) +
  geom_bar(stat = "identity", fill = "firebrick") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "P-valores por predictor",
    x = "Predictor",
    y = "P-valor"
  ) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") + # Línea de significancia
  geom_text(aes(label = round(`p-valor`, 3)), hjust = -0.1, size = 4) +
  ylim(0, max(tabla_publicacion$`p-valor`, na.rm = TRUE) * 1.2)


coef_modelo <- coef(modelo_final)  
coef_mean <- rowMeans(coef_modelo)
coef_se <- apply(coef_modelo, 1, sd)


coef_df <- data.frame(
  Predictora = names(coef_mean),
  Coeficiente = coef_mean,
  Inferior = coef_mean - 1.96 * coef_se,
  Superior = coef_mean + 1.96 * coef_se
)


coef_df <- coef_df[-6,]
coef_df <- coef_df[-1,] 


coef_df$Predictora <- gsub("\\(Intercept\\)", "Intercepto", coef_df$Predictora)

ggplot(coef_df, aes(x = reorder(Predictora, Coeficiente), y = Coeficiente, fill = Coeficiente > 0)) +
  geom_col(show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("firebrick3", "forestgreen")) +
  labs(
    title = "Dirección del efecto de cada predictor",
    x = "Predictor",
    y = "Coeficiente estimado"
  )







##### graficar los residuos de las variables para ver las tendencias


library(patchwork)
library(viridis)

residuales <- modelo_final$pgls.residuals

pca_resid <- prcomp(residuales, scale. = TRUE)
Y1_resid <- pca_resid$x[, 1]

df_residuos <- as.data.frame(modelo_final$data)
df_residuos$Morpho_PC1_resid <- Y1_resid

variables_modelo <- c("foraging", "bmass", "centroid_size", "range", "maxelev", "abs_centlatitude")

paleta_extendida <- c(paleta_foraging, "orange", "purple", "brown", "grey40")

plots_lista_resid <- list()

for (var in variables_modelo) {
  
  df_temp <- df_residuos
  
  if (var == "range") {
    limite_range <- quantile(df_temp$range, 0.95, na.rm = TRUE)
    df_temp <- df_temp %>% filter(range <= limite_range)
  }
  
  if (var == "centroid_size") {
    limite_size <- quantile(df_temp$centroid_size, 0.99, na.rm = TRUE)
    df_temp <- df_temp %>% filter(centroid_size <= limite_size)
  }
  
  if (is.factor(df_temp[[var]]) || is.character(df_temp[[var]])) {

    color_palette <- if (var == "foraging") {
      paleta_foraging
    } else {
      paleta_extendida
    }
    
    p <- ggplot(df_temp, aes_string(x = var, y = "Morpho_PC1_resid", fill = var, color = var)) +
      geom_boxplot(alpha = 0.6, outlier.shape = NA) +
      geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.6) +
      scale_fill_manual(values = color_palette) +
      scale_color_manual(values = color_palette) +
      theme_minimal(base_size = 14) +
      labs(
        title = paste("Residual morphology vs", var),
        x = var,
        y = "Morphological residuals (PC1)"
      )
    
  } else {
  
    p <- ggplot(df_temp, aes_string(x = var, y = "Morpho_PC1_resid")) +
      geom_point(size = 2, alpha = 0.6, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "purple", size = 1) +
      theme_minimal(base_size = 14) +
      labs(
        title = paste("Residual morphology vs", var),
        x = var,
        y = "Morphological residuals (PC1)"
      )
  }
  

  plots_lista_resid[[var]] <- p
}

wrap_plots(plots_lista_resid)






###### Grafico 2, datos crudos

predicciones <- modelo_final$pgls.fitted

pca_pred <- prcomp(predicciones, scale. = TRUE)
Y1 <- pca_pred$x[, 1]

df_predictores <- as.data.frame(modelo_final$data)
df_predictores$Morpho_PC1 <- Y1

variables_modelo <- c("foraging", "bmass", "centroid_size", "range", "maxelev", "abs_centlatitude")

paleta_extendida <- c(paleta_foraging, "orange", "purple", "brown", "grey40")

plots_lista <- list()

for (var in variables_modelo) {
  
  df_temp <- df_predictores
  
  if (var == "range") {
    limite_range <- quantile(df_temp$range, 0.95, na.rm = TRUE)
    df_temp <- df_temp %>% filter(range <= limite_range)
  }
  
  if (var == "centroid_size") {
    limite_size <- quantile(df_temp$centroid_size, 0.99, na.rm = TRUE)
    df_temp <- df_temp %>% filter(centroid_size <= limite_size)
  }
  
  if (is.factor(df_temp[[var]]) || is.character(df_temp[[var]])) {

    color_palette <- if (var == "foraging") paleta_foraging else paleta_extendida
    
    p <- ggplot(df_temp, aes_string(x = var, y = "Morpho_PC1", fill = var, color = var)) +
      geom_boxplot(alpha = 0.6, outlier.shape = NA) +
      geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.6) +
      scale_fill_manual(values = color_palette) +
      scale_color_manual(values = color_palette) +
      theme_light(base_size = 14) +
      theme(panel.grid = element_blank()) +
      labs(
        title = paste(var),
        x = var,
        y = "Humeral shape (PC1)"
      )
    
  } else {

    p <- ggplot(df_temp, aes_string(x = var, y = "Morpho_PC1")) +
      geom_point(size = 2, alpha = 0.6, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "purple", size = 1) +
      theme_light(base_size = 14) +
      theme(panel.grid = element_blank()) +
      labs(
        title = paste(var),
        x = var,
        y = "Humeral shape (PC1)"
      )
  }
  
  plots_lista[[var]] <- p
}

wrap_plots(plots_lista)





### Grafico 3, sin outliers

predicciones <- modelo_final$pgls.fitted

pca_pred <- prcomp(predicciones, scale. = TRUE)
Y1 <- pca_pred$x[, 1]

df_predictores <- as.data.frame(modelo_final$data)
df_predictores$Morpho_PC1 <- Y1

variables_modelo <- c("foraging", "bmass", "centroid_size", "range", "maxelev", "abs_centlatitude")

paleta_extendida <- c(paleta_foraging, "orange", "purple", "brown", "grey40")

plots_lista <- list()

for (var in variables_modelo) {

  df_temp <- df_predictores

  if (var == "range") {
    limite_range <- quantile(df_temp$range, 0.95, na.rm = TRUE)
    df_temp <- df_temp %>% filter(range <= limite_range)
  }

  if (var == "centroid_size") {
    limite_size <- quantile(df_temp$centroid_size, 0.99, na.rm = TRUE)
    df_temp <- df_temp %>% filter(centroid_size <= limite_size)
  }

  if (var == "abs_centlatitude") {
    df_temp <- df_temp %>% filter(abs_centlatitude <= 30)
  }
  
  if (is.factor(df_temp[[var]]) || is.character(df_temp[[var]])) {
    color_palette <- if (var == "foraging") paleta_foraging else paleta_extendida
    
    p <- ggplot(df_temp, aes_string(x = var, y = "Morpho_PC1", fill = var, color = var)) +
      geom_boxplot(alpha = 0.6, outlier.shape = NA) +
      geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.6) +
      scale_fill_manual(values = color_palette) +
      scale_color_manual(values = color_palette) +
      theme_light(base_size = 14) +
      theme(panel.grid = element_blank()) +
      labs(
        title = paste(var),
        x = var,
        y = "Humeral shape (PC1)"
      )
    
  } else {

    p <- ggplot(df_temp, aes_string(x = var, y = "Morpho_PC1")) +
      geom_point(size = 2, alpha = 0.6, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "purple", size = 1) +
      theme_light(base_size = 14) +
      theme(panel.grid = element_blank()) +
      labs(
        title = paste(var),
        x = var,
        y = "Humeral shape (PC1)"
      )
  }
  
  plots_lista[[var]] <- p
}

wrap_plots(plots_lista)


dev.off()


#PCA - differences

# Phylomorphospace con colores por grupo de forrajeo

plot.res <- gm.prcomp(shape_clean, phy = tree_clean)

ind_names <- dimnames(shape_clean)[[3]]

foraging_plot <- metadata_clean[ind_names, "foraging"]
phylo.par <- "tip.labels"
phylo.par <- as.list(phylo.par)

plot(plot.res,
     phylo = TRUE,
     pch = 21,
     bg = paleta_foraging[foraging_plot],
     cex = 1.5,
     phylo.par = list(
       tip.labels = FALSE,
       node.labels = FALSE,
       node.cex = 0.1
     ))



#Phylogenetic signal

PS.shape <- physignal(gdf$shape, gdf$tree, print.progress = FALSE)
summary(PS.shape)
plot(PS.shape)

K_obs <- mean(PS.shape$K)

K_null_means <- rowMeans(matrix(PS.shape$random.K, ncol = length(PS.shape$K), byrow = TRUE))

Z_effect <- (K_obs - mean(K_null_means)) / sd(K_null_means)

cat("📏 Tamaño del efecto Z (global):", round(Z_effect, 3), "\n")


# Integración morfológica con corrección filogenética

land.gps <- c("A","A","A","A","A","A","A","A","A","B","B","B","B","B","B","B","B","B","B","B","B"
              ,"A","A","A","A")

PLS.Y <- phylo.integration(A = gdf$shape, partition.gp = land.gps, phy = tree_clean, print.progress = FALSE)
summary(PLS.Y)
plot(PLS.Y)


pls_df <- data.frame(
  Block1 = scores_block1[, 1],
  Block2 = scores_block2[, 1]
)


pls_df_clean <- pls_df %>%
  mutate(
    z_Block1 = scale(Block1),
    z_Block2 = scale(Block2)
  ) %>%
  filter(abs(z_Block1) < 3, abs(z_Block2) < 3)

ggplot(pls_df_clean, aes(x = Block1, y = Block2)) +
  geom_point(size = 2, alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "purple", size = 1) +
  labs(
    title = "Phenotypic Integration of Humeral Regions",
    x = "PLS Axis 1 (Proximal region)",
    y = "PLS Axis 1 (Distal region)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )


## Tasas evolutivas


names(foraging) <- rownames(metadata_clean)
names(habitat) <- rownames(metadata_clean)

# Modelo 1: tasas por foraging
ER_foraging <- compare.evol.rates(
  A = shape_clean,
  phy = tree_clean,
  gp = foraging,
  iter = 999,
  method = "permutation",
  print.progress = FALSE
)

summary(ER_foraging)
plot(ER_foraging, main = "Tasa de evolución por estrategia de forrajeo")

# variables continuas (ecológicas)
dimnames(shape_clean)[[3]] <- rownames(metadata_clean)

species_list <- tree_clean$tip.label


species_factor <- factor(dimnames(shape_clean)[[3]])

library(geomorph)
species_means <- arrayspecs(
  aggregate(two.d.array(shape_clean), by = list(Species = species_factor), FUN = mean)[, -1],
  p = dim(shape_clean)[1], 
  k = dim(shape_clean)[2]
)
dimnames(species_means)[[3]] <- levels(species_factor)

species_matrix <- two.d.array(species_means[, , species_list])


species_means <- shape_clean[, , tree_clean$tip.label]

species_matrix <- two.d.array(species_means)
rownames(species_matrix) <- tree_clean$tip.label


anc_coords <- apply(species_matrix, 2, function(col) {
  names(col) <- rownames(species_matrix)
  fastAnc(tree_clean, col)
})


anc_coords <- as.matrix(anc_coords)

ntips <- length(tree_clean$tip.label)
nnodes <- tree_clean$Nnode

rownames(species_matrix) <- 1:ntips
rownames(anc_coords) <- (ntips + 1):(ntips + nnodes)


all_shapes <- rbind(species_matrix, anc_coords)


branch_change <- sapply(1:nrow(tree_clean$edge), function(i) {
  from <- tree_clean$edge[i, 1]
  to <- tree_clean$edge[i, 2]
  dist(rbind(all_shapes[as.character(from), ], all_shapes[as.character(to), ]))
})



tip_edges_idx <- which(tree_clean$edge[, 2] <= length(tree_clean$tip.label))
tip_nodes <- tree_clean$edge[tip_edges_idx, 2]
tip_labels <- tree_clean$tip.label[tip_nodes]

branching_times <- branching.times(tree_clean)
total_time <- max(branching_times)
tip_ages <- node.depth.edgelength(tree_clean)[tip_nodes]
divergence_time <- total_time - tip_ages

morph_change <- branch_change[tip_edges_idx]
morph_rate <- morph_change / divergence_time

valid_idx <- which(is.finite(morph_rate) & morph_rate > 0 & divergence_time > 0)


plot_df <- data.frame(
  Species = tip_labels[valid_idx],
  MorphRate = morph_rate[valid_idx],
  bmass = metadata_clean[tip_labels[valid_idx], "bmass"],
  range = metadata_clean[tip_labels[valid_idx], "Range.Size"],
  centroid_size = metadata_clean[tip_labels[valid_idx], "centroid_size"],
  maxelev = metadata_clean[tip_labels[valid_idx], "maxelev"],
  abs_centlatitude = metadata_clean[tip_labels[valid_idx], "abs_centlatitude"]
)


q975_rate <- quantile(plot_df$MorphRate, 0.97, na.rm = TRUE)
plot_df_clean <- plot_df %>%
  filter(MorphRate <= q975_rate)


predictors <- c("bmass", "abs_centlatitude", "range", "centroid_size", "maxelev")

modelos_continuos <- list()

tabla_resultados <- data.frame(
  Predictor = character(),
  Estimate = numeric(),
  StdError = numeric(),
  t_value = numeric(),
  p_value = numeric(),
  R_squared = numeric(),
  stringsAsFactors = FALSE
)

library(mgcv)
# Ajuste del modelo GAM
modelo_gam <- gam(
  MorphRate ~ s(abs_centlatitude) + s(centroid_size) + s(bmass) + s(range) + s(maxelev),
  data = plot_df_clean,
  method = "REML"
)

summary(modelo_gam)

plot(modelo_gam, pages = 1, residuals = TRUE, pch = 16, cex = 0.8, rug = TRUE)
gam.check(modelo_gam)





##### Graficar tasas evolutivas


plot_df_clean <- plot_df_clean %>% filter(MorphRate <= quantile(MorphRate, 0.975, na.rm = TRUE))


plot_df_long <- plot_df_clean %>%
  pivot_longer(cols = c(bmass, range, centroid_size, maxelev, abs_centlatitude),
               names_to = "Variable",
               values_to = "Valor")

plot_continuas <- ggplot(plot_df_long, aes(x = Valor, y = MorphRate)) +
  geom_point(color = "steelblue", size = 1, alpha = 0.6) +
  geom_smooth(method = "gam", se = TRUE, color = "purple", size = 1) +
  facet_wrap(~ Variable, scales = "free_x", ncol = 2) +
  theme_light(base_size = 14) +
  theme(panel.grid = element_blank()) +
  labs(
    title = "Tasa de evolución morfológica vs. predictores continuos",
    x = "Predictor value",
    y = "Humeral morphology evolutionary rate"
  )

print(plot_continuas)




###### Graficar tasas en filogenia

branch_rates <- branch_change / tree_clean$edge.length
branch_rates <- rank(branch_rates) / length(branch_rates)


colors <- viridis(100)
color_index <- as.numeric(cut(branch_rates, breaks = 100))
edge_colors <- colors[color_index]


layout(matrix(c(1, 2), 1, 2), widths = c(4, 1))


par(mar = c(1, 1, 2, 1), oma = c(0, 0, 4, 0))


plot(tree_clean,
     type = "fan",
     edge.color = edge_colors,
     edge.width = 2.5,
     cex = 0.6,
     label.offset = 0.01,
     no.margin = TRUE)


mtext("Tasa de evolución morfológica por rama", side = 3, outer = TRUE, line = 1.5, cex = 1.2)

par(new = TRUE, mar = c(5, 1, 5, 4))
image.plot(legend.only = TRUE,
           zlim = range(branch_rates),
           col = colors,
           legend.shrink = 0.5,
           legend.width = 1.2,
           smallplot = c(0.87, 0.89, 0.2, 0.8),
           legend.args = list(text = "Tasa morfológica", side = 3, line = 1))






dev.off()


## Grafica 2

add_alpha <- function(col, alpha = 1) {
  rgb_col <- col2rgb(col, alpha = TRUE)
  rgb(rgb_col[1, ], rgb_col[2, ], rgb_col[3, ], alpha = alpha * 255, maxColorValue = 255)
}



tree_scaled <- tree_clean
tree_scaled$edge.length <- tree_scaled$edge.length * 0.55

branch_rates <- branch_change / tree_clean$edge.length
branch_rates <- rank(branch_rates) / length(branch_rates)
colors <- viridis(100)
color_index <- as.numeric(cut(branch_rates, breaks = 100))
edge_colors <- colors[color_index]

tip_labels <- tree_scaled$tip.label
species_index <- match(tip_labels, rownames(metadata_clean))
tip_groups_foraging <- metadata_clean$foraging[species_index]
tip_colors_foraging <- paleta_foraging[tip_groups_foraging]

scale_log <- function(var) {
  log_var <- log1p(var)
  scaled <- (log_var - min(log_var, na.rm = TRUE)) / (max(log_var, na.rm = TRUE) - min(log_var, na.rm = TRUE))
  return(0.8 + scaled * 3.5)
}

elev_lengths <- scale_log(metadata_clean$maxelev[species_index])
lat_lengths <- scale_log(abs(metadata_clean$`Centroid.Latitude`[species_index]))
range_lengths <- scale_log(metadata_clean$Range.Size[species_index])
bmass_lengths <- scale_log(metadata_clean$bmass[species_index])

par(mar = c(1, 1, 1, 1), oma = c(0, 0, 4, 0))
tree_radius <- max(node.depth.edgelength(tree_scaled)) * 1.8

plot(tree_scaled,
     type = "fan",
     edge.color = edge_colors,
     edge.width = 2.5,
     cex = 0.6,
     label.offset = 0.015,
     no.margin = TRUE,
     x.lim = c(-tree_radius, tree_radius),
     y.lim = c(-tree_radius, tree_radius))

tree_layout <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_coords <- cbind(tree_layout$xx[1:Ntip(tree_scaled)],
                    tree_layout$yy[1:Ntip(tree_scaled)])
tip_vectors <- tip_coords / sqrt(rowSums(tip_coords^2))

angles <- atan2(tip_coords[, 2], tip_coords[, 1])
tip_order <- order(angles)

tip_coords <- tip_coords[tip_order, ]
tip_vectors <- tip_vectors[tip_order, ]
tip_colors_foraging <- tip_colors_foraging[tip_order]
elev_lengths <- elev_lengths[tip_order]
lat_lengths <- lat_lengths[tip_order]
range_lengths <- range_lengths[tip_order]
bmass_lengths <- bmass_lengths[tip_order]


expanded_coords <- tip_coords + tip_vectors * 4
points(expanded_coords, pch = 21, bg = tip_colors_foraging, col = "black", cex = 1.5, lwd = 0.6)


tangents <- cbind(-tip_vectors[, 2], tip_vectors[, 1])
distance_from_center <- 17
bar_center <- tip_vectors * distance_from_center

offsets <- c(-0.45, -0.15, 0.15, 0.45)
bar_lengths_raw <- list(elev_lengths, lat_lengths, range_lengths, bmass_lengths)
bar_colors <- c("gray30", "dodgerblue3", "deeppink3", "goldenrod3")

for (i in seq_along(offsets)) {
  segments(
    bar_center[,1] + tangents[,1] * offsets[i],
    bar_center[,2] + tangents[,2] * offsets[i],
    bar_center[,1] + tangents[,1] * offsets[i] + tip_vectors[,1] * bar_lengths_raw[[i]],
    bar_center[,2] + tangents[,2] * offsets[i] + tip_vectors[,2] * bar_lengths_raw[[i]],
    col = bar_colors[i],
    lwd = 2.2
  )
}


mtext("Tasa de evolución morfológica por rama", side = 3, outer = TRUE, line = 1.5, cex = 1.2)

par(new = TRUE, mar = c(5, 1, 5, 4))
image.plot(legend.only = TRUE,
           zlim = range(branch_rates),
           col = colors,
           legend.shrink = 0.5,
           legend.width = 1.2,
           smallplot = c(0.87, 0.89, 0.2, 0.8),
           legend.args = list(text = "Tasa morfológica", side = 3, line = 1))

legend("topleft",
       legend = names(paleta_foraging),
       pt.bg = paleta_foraging,
       pch = 21,
       bty = "n",
       title = "Forrajeo",
       cex = 0.9)

legend(x = par("usr")[1] + 0.001 * diff(par("usr")[1:2]),
       y = par("usr")[4] - 0.05 * diff(par("usr")[3:4]),
       legend = c("Elevación", "Latitud", "Tamaño del rango", "Masa corporal (bmass)"),
       col = bar_colors,
       lwd = 2.2,
       title = "Variables continuas",
       bty = "n",
       cex = 0.9)





dev.off()




#### figura de frecuencias de tasas evolutivas

log_branch_rates <- log(branch_rates)


branch_bins <- cut(log_branch_rates, breaks = 30) 
palette_colors <- viridis(length(levels(branch_bins)), option = "D")
names(palette_colors) <- levels(branch_bins)

df_hist <- data.frame(
  log_rate = log_branch_rates,
  bin = branch_bins
)

ggplot(df_hist, aes(x = log_rate, fill = bin)) +
  geom_histogram(bins = 30, color = NA) +
  scale_fill_manual(values = palette_colors, guide = "none") +
  theme_minimal(base_size = 14) +
  labs(x = "log(Mean Rate)", y = NULL)







####################


# PCA
pca_res <- gm.prcomp(shape_clean)
coords_pca <- pca_res$x[, 1:2]
specimens <- dimnames(shape_clean)[[3]]
rownames(coords_pca) <- specimens

group_raw <- metadata_clean[specimens, "foraging"]
group_factor <- factor(group_raw, levels = c("ter", "tra", "opp"))

tip_colors <- paleta_foraging[group_factor]

# Gráfico con phylomorphospace
phylomorphospace(tree_clean, coords_pca,
                 label = "off",
                 node.size = c(0, 0.5),
                 xlab = "PC1",
                 ylab = "PC2",
                 main = "Espacio morfológico con estructura filogenética",
                 control = list(col = tip_colors))

tiplabels(pch = 21, bg = tip_colors, cex = 1.2)

legend("topleft",
       legend = levels(group_factor),
       pch = 21,
       pt.bg = paleta_foraging,
       title = "Forrajeo",
       bty = "n")


##### dotplot

library(gridExtra)

foraging_df <- data.frame(
  Grupo = c("ter", "opp", "tra"),
  Tasa = c(1.541896e-05, 1.716319e-05, 3.876088e-05)
)
foraging_df$Grupo <- factor(foraging_df$Grupo, levels = c("ter", "opp", "tra"))

plot_dot <- ggplot(foraging_df, aes(x = Grupo, y = Tasa, color = Grupo)) +
  geom_point(size = 5) +
  geom_segment(aes(xend = Grupo, y = 0, yend = Tasa), size = 1.2) +
  scale_color_manual(values = paleta_foraging) +
  theme_minimal() +
  labs(title = "Tasa evolutiva por grupo de forrajeo", x = "Foraging", y = "Morphological evolutionary rate") +
  theme(legend.position = "none")


groups <- factor(metadata_clean$foraging)
group_levels <- levels(groups)


length(group_levels)        
dim(shape_clean)[3]         
length(groups)              

names(groups) <- dimnames(shape_clean)[[3]]

groups <- groups[dimnames(shape_clean)[[3]]]

group_levels <- levels(groups)

consensus_shapes <- lapply(group_levels, function(gr) {
  gp_coords <- shape_clean[, , groups == gr]
  mshape(gp_coords)
})


outline_df <- do.call(rbind, lapply(1:length(group_levels), function(i) {
  coords <- consensus_shapes[[i]]
  data.frame(
    X = coords[, 1],
    Y = coords[, 2],
    Grupo = group_levels[i],
    Landmark = 1:nrow(coords)
  )
}))

plot_shapes <- ggplot(outline_df, aes(x = X, y = Y, group = Grupo, color = Grupo)) +
  geom_path(size = 1.2) +
  coord_equal() +
  theme_void() +
  scale_color_manual(values = paleta_foraging) +
  facet_wrap(~ Grupo) +
  labs(title = "Forma consenso por grupo")

grid.arrange(plot_dot, plot_shapes, ncol = 2, widths = c(1, 2))

outline_df <- do.call(rbind, lapply(1:length(group_levels), function(i) {
  coords <- consensus_shapes[[i]]
  data.frame(
    X = coords[, 1],
    Y = coords[, 2],
    Grupo = group_levels[i],
    Landmark = 1:nrow(coords)
  )
}))


ggplot(outline_df, aes(x = X, y = Y, group = Grupo, color = Grupo)) +
  geom_path(size = 1.2) +
  coord_equal() +
  theme_minimal() +
  scale_color_manual(values = paleta_foraging) +
  labs(title = "Formas consenso por grupo (superpuestas)",
       x = NULL, y = NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )






##### CVA corregido por centroid size
library(Morpho)

calc_centroid_size <- function(coords) {
  apply(coords, 3, function(x) sqrt(sum((x - colMeans(x))^2)))
}

# Calcular centroid size para cada espécimen
cs <- calc_centroid_size(shape_clean)

log_shape <- two.d.array(shape_clean)

shape_resid_model <- procD.lm(log_shape ~ log(cs), iter = 999, print.progress = FALSE)

shape_resid <- arrayspecs(shape_resid_model$residuals, 
                          p = dim(shape_clean)[1], 
                          k = dim(shape_clean)[2])

foraging_group <- factor(metadata_clean$foraging)

cva_res <- CVA(shape_resid, 
               group = foraging_group, 
               cv = FALSE, 
               plot = TRUE)

# permutaciones
fit_size <- procD.lm(shape_clean ~ cs)

shape_cva <- arrayspecs(residuals(fit_size), 
                        p = dim(shape_clean)[1], 
                        k = dim(shape_clean)[2])



set.seed(123)
n_perm <- 999
observed_acc <- mean(cva_res$class == cva_res$groups)

perm_accs <- numeric(n_perm)

for (i in 1:n_perm) {
  perm_groups <- sample(cva_res$groups)
  perm_cva <- CVA(shape_cva, groups = perm_groups, lda = TRUE)
  
  perm_accs[i] <- mean(perm_cva$class == perm_groups)
}

p_val <- mean(perm_accs >= observed_acc)

cat("🔁 Permutation test (", n_perm, " perm):\n", sep = "")
cat("Observed classification accuracy:", round(observed_acc, 4), "\n")
cat("Empirical p-value:", round(p_val, 4), "\n")


scores <- as.data.frame(cva_res$CVscores)
colnames(scores)[1:2] <- c("CV1", "CV2")
scores$Group <- cva_res$groups

ggplot(scores, aes(x = CV1, y = CV2, color = Group, fill = Group)) +
  geom_point(size = 3, shape = 21, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.9, linetype = "dashed", size = 1) +
  scale_color_manual(values = paleta_foraging) +
  scale_fill_manual(values = paleta_foraging) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Espacio Canónico (CVA) según tipo de forrajeo",
    x = "CV1", y = "CV2",
    color = "Forrajeo", fill = "Forrajeo"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )


dim(shape_clean)
length(foraging)

foraging_factor <- as.factor(foraging)
names(foraging_factor) <- dimnames(shape_clean)[[3]]

shape_2d <- two.d.array(shape_clean)

cva_res <- CVA(shape_2d, groups = foraging_factor, lda = TRUE)

if (nrow(cva_res$CVscores) == 0) stop("❌ CVA no generó scores. Revisa los datos.")

scores <- data.frame(
  CV1 = cva_res$CVscores[, 1],
  CV2 = cva_res$CVscores[, 2],
  Group = foraging_factor
)

centroids <- aggregate(cbind(CV1, CV2) ~ Group, data = scores, FUN = mean)

ggplot(scores, aes(x = CV1, y = CV2, color = Group, fill = Group)) +
  geom_point(size = 3, shape = 21, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.9, linetype = "dashed", size = 1) +
  geom_point(data = centroids, aes(x = CV1, y = CV2),
             shape = 23, size = 4, fill = "white", color = "black", stroke = 1.2) +
  scale_color_manual(values = paleta_foraging) +
  scale_fill_manual(values = paleta_foraging) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Espacio Canónico (CVA) según tipo de forrajeo",
    x = "CV1", y = "CV2",
    color = "Forrajeo", fill = "Forrajeo"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

### grafica 2  CVA

scores <- cva_res$CVscores
rownames(scores) <- names(foraging_factor)

common_tips <- intersect(rownames(scores), tree_clean$tip.label)
tree_cv <- drop.tip(tree_clean, setdiff(tree_clean$tip.label, common_tips))
scores_cv <- scores[tree_cv$tip.label, 1:2]
groups_cv <- foraging_factor[tree_cv$tip.label]

groups_cv <- factor(groups_cv, levels = names(paleta_foraging))

group_colors <- setNames(paleta_foraging, names(paleta_foraging))
tip_colors <- group_colors[as.character(groups_cv)]

phylomorphospace(tree_cv, scores_cv,
                 label = "off", node.size = c(0, 0),
                 xlab = "CV1", ylab = "CV2",
                 main = "Filogenia sobre espacio canónico (CVA)")

points(scores_cv, pch = 21, bg = tip_colors, col = "black", cex = 1.8, lwd = 0.7)

library(ellipse)
library(MASS)

for (grp in levels(groups_cv)) {
  grp_scores <- scores_cv[groups_cv == grp, , drop = FALSE]
  
  if (nrow(grp_scores) > 2) {
    cov_grp <- cov.trob(grp_scores)$cov
    center <- colMeans(grp_scores)
    elip_coords <- ellipse::ellipse(x = cov_grp, centre = center, level = 0.90)
    
    lines(elip_coords, col = group_colors[grp], lty = 2, lwd = 2)
  }
}

legend("topright",
       legend = names(paleta_foraging),
       pt.bg = paleta_foraging,
       pch = 21,
       title = "Forrajeo",
       bty = "n",
       cex = 0.9)