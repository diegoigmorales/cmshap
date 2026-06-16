# =========================================================
# 05_boruta.R
# Relevancia preliminar de atributos con Boruta en Swissmetro
# =========================================================

# -----------------------------
# 0) Paquetes
# -----------------------------
library(Boruta)

set.seed(2026)

# -----------------------------
# 1) Rutas
# -----------------------------
DATA_FILE <- "data/processed/swissmetro_clean.csv"
OUTPUT_DIR <- "output"
TABLE_DIR <- file.path(OUTPUT_DIR, "tablas")
FIGURE_DIR <- file.path(OUTPUT_DIR, "figuras")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2) Cargar y preparar datos
# -----------------------------
sm <- read.csv(DATA_FILE)

sm$choice <- factor(
  sm$choice,
  levels = c(1, 2, 3),
  labels = c("train", "swissmetro", "car")
)

factor_vars <- c(
  "age", "ga", "luggage", "male", "income",
  "purpose", "first", "group", "sm_seats"
)
sm[factor_vars] <- lapply(sm[factor_vars], factor)

boruta_vars <- c(
  "choice",
  "train_tt", "train_co", "train_he",
  "sm_tt", "sm_co", "sm_he", "sm_seats",
  "car_tt", "car_co",
  "age", "ga", "luggage",
  "male", "income", "purpose", "first", "group"
)

boruta_data <- na.omit(sm[boruta_vars])

# -----------------------------
# 3) Ejecutar Boruta
# -----------------------------
boruta_fit <- Boruta(
  choice ~ .,
  data = boruta_data,
  maxRuns = 200,
  doTrace = 1
)

boruta_final <- TentativeRoughFix(boruta_fit)

# -----------------------------
# 4) Guardar tabla de importancia
# -----------------------------
grupo_mnl <- c(
  train_tt = "time", sm_tt = "time", car_tt = "time",
  train_co = "cost", sm_co = "cost", car_co = "cost",
  train_he = "headway", sm_he = "headway",
  sm_seats = "seats",
  age = "age", ga = "ga", luggage = "luggage",
  male = "socioeconomic", income = "socioeconomic",
  purpose = "trip_context", first = "trip_context", group = "trip_context"
)

boruta_importance <- attStats(boruta_final)
boruta_importance <- data.frame(
  variable = rownames(boruta_importance),
  boruta_importance,
  row.names = NULL
)
boruta_importance$grupo_mnl <- unname(grupo_mnl[boruta_importance$variable])
boruta_importance$grupo_mnl[is.na(boruta_importance$grupo_mnl)] <- "other"

orden <- order(
  boruta_importance$decision != "Confirmed",
  -boruta_importance$medianImp
)
boruta_importance <- boruta_importance[orden, ]

write.csv(
  boruta_importance,
  file.path(TABLE_DIR, "boruta_importancia_swissmetro.csv"),
  row.names = FALSE
)

# -----------------------------
# 5) Guardar grafico
# -----------------------------
png(
  filename = file.path(FIGURE_DIR, "boruta_importancia_swissmetro.png"),
  width = 1300,
  height = 800,
  res = 120
)
plot(
  boruta_final,
  sort = TRUE,
  las = 2,
  cex.axis = 0.75,
  main = "Boruta: importancia preliminar de atributos"
)
dev.off()

cat("=====================================================\n")
cat("Boruta terminado\n")
cat("=====================================================\n")
cat("Observaciones usadas:", nrow(boruta_data), "\n\n")

cat("Variables confirmadas:\n")
print(boruta_importance[boruta_importance$decision == "Confirmed", ])
cat("\n")

cat("Tabla guardada en:", file.path(TABLE_DIR, "boruta_importancia_swissmetro.csv"), "\n")
cat("Grafico guardado en:", file.path(FIGURE_DIR, "boruta_importancia_swissmetro.png"), "\n")
