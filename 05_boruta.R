# =========================================================
# 05_boruta.R
# Relevancia preliminar de atributos con Boruta en Swissmetro
# =========================================================

# 1. Cargar paquete
library(Boruta)

set.seed(2026)

# 2. Definir rutas
archivo_entrada <- "swissmetro_limpio.csv"
directorio_salida <- "output"
directorio_tablas <- file.path(directorio_salida, "tablas")
directorio_figuras <- file.path(directorio_salida, "figuras")

dir.create(directorio_tablas, recursive = TRUE, showWarnings = FALSE)
dir.create(directorio_figuras, recursive = TRUE, showWarnings = FALSE)

# 3. Cargar datos resultantes de 02_depurar.R
if (!file.exists(archivo_entrada)) {
  source("02_depurar.R")
}

df <- read.csv(archivo_entrada, stringsAsFactors = FALSE)

# 4. Preparar variables para Boruta
df$CHOICE <- factor(
  df$CHOICE,
  levels = c(1, 2, 3),
  labels = c("train", "swissmetro", "car")
)

variables_factor <- c("AGE", "GA", "SM_SEATS", "LUGGAGE")
df[variables_factor] <- lapply(df[variables_factor], factor)

variables_boruta <- c(
  "CHOICE",
  "CAR_TT", "TRAIN_TT", "SM_TT",
  "CAR_CO", "TRAIN_CO", "SM_CO",
  "AGE", "GA", "SM_SEATS", "LUGGAGE"
)

datos_boruta <- na.omit(df[variables_boruta])

# 5. Ejecutar Boruta
ajuste_boruta <- Boruta(
  CHOICE ~ .,
  data = datos_boruta,
  maxRuns = 200,
  doTrace = 1
)

ajuste_boruta_final <- TentativeRoughFix(ajuste_boruta)

# 6. Guardar tabla de importancia
grupo_mnl <- c(
  CAR_TT = "TT",
  TRAIN_TT = "TT",
  SM_TT = "TT",
  CAR_CO = "COST",
  TRAIN_CO = "COST",
  SM_CO = "COST",
  AGE = "AGE",
  GA = "GA",
  SM_SEATS = "SEATS",
  LUGGAGE = "LUGGAGE"
)

importancia_boruta <- attStats(ajuste_boruta_final)
importancia_boruta <- data.frame(
  variable = rownames(importancia_boruta),
  importancia_boruta,
  row.names = NULL
)

importancia_boruta$grupo_mnl <- unname(grupo_mnl[importancia_boruta$variable])
importancia_boruta$grupo_mnl[is.na(importancia_boruta$grupo_mnl)] <- "other"

orden <- order(
  importancia_boruta$decision != "Confirmed",
  -importancia_boruta$medianImp
)

importancia_boruta <- importancia_boruta[orden, ]

write.csv(
  importancia_boruta,
  file.path(directorio_tablas, "boruta_importancia_swissmetro.csv"),
  row.names = FALSE
)

# 7. Guardar grafico
png(
  filename = file.path(directorio_figuras, "boruta_importancia_swissmetro.png"),
  width = 1300,
  height = 800,
  res = 120
)

plot(
  ajuste_boruta_final,
  sort = TRUE,
  las = 2,
  cex.axis = 0.75,
  main = "Boruta: importancia preliminar de atributos"
)

dev.off()

# 8. Imprimir resultados principales
cat("=====================================================\n")
cat("Boruta terminado\n")
cat("=====================================================\n")
cat("Observaciones usadas:", nrow(datos_boruta), "\n\n")

cat("Variables confirmadas:\n")
print(importancia_boruta[importancia_boruta$decision == "Confirmed", ])
cat("\n")

cat("Tabla guardada en:", file.path(directorio_tablas, "boruta_importancia_swissmetro.csv"), "\n")
cat("Grafico guardado en:", file.path(directorio_figuras, "boruta_importancia_swissmetro.png"), "\n")
