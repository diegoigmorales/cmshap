# =========================================================
# 02_depurar_swissmetro.R
# Depuración de Swissmetro para la tesis
# =========================================================

# -----------------------------
# 0) Paquetes
# -----------------------------
library(readr)
library(dplyr)
library(janitor)

# -----------------------------
# 1) Rutas
# -----------------------------
archivo_entrada <- "swissmetro.dat"
dir.create("data", showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

archivo_salida_csv <- "data/processed/swissmetro_clean.csv"
archivo_salida_rds <- "data/processed/swissmetro_clean.rds"

# -----------------------------
# 2) Cargar datos
# -----------------------------
# El archivo viene tabulado
sm_raw <- read_tsv(archivo_entrada, show_col_types = FALSE) |>
  clean_names()

cat("=====================================================\n")
cat("Swissmetro: auditoría inicial\n")
cat("=====================================================\n")
cat("Observaciones originales:", nrow(sm_raw), "\n")
cat("Número de columnas:", ncol(sm_raw), "\n\n")

cat("Nombres de variables:\n")
print(names(sm_raw))
cat("\n")

cat("Valores perdidos por variable:\n")
print(colSums(is.na(sm_raw)))
cat("\n")

cat("Duplicados exactos:", sum(duplicated(sm_raw)), "\n\n")

# -----------------------------
# 3) Verificación rápida de variables clave
# -----------------------------
variables_clave <- c(
  "id", "choice",
  "train_av", "car_av", "sm_av",
  "train_tt", "train_co",
  "sm_tt", "sm_co",
  "car_tt", "car_co",
  "age", "ga", "sm_seats", "luggage"
)

faltan <- setdiff(variables_clave, names(sm_raw))

if (length(faltan) > 0) {
  stop(
    paste(
      "Faltan estas variables clave en la base:",
      paste(faltan, collapse = ", ")
    )
  )
}

# -----------------------------
# 4) Depuración principal
# -----------------------------
# Criterio:
# - eliminar casos sin elección válida
# - dejar solo observaciones con las 3 alternativas disponibles
sm_clean <- sm_raw |>
  filter(
    !is.na(choice),
    choice != 0,
    train_av == 1,
    car_av == 1,
    sm_av == 1
  )

cat("=====================================================\n")
cat("Después del filtro principal\n")
cat("=====================================================\n")
cat("Observaciones restantes:", nrow(sm_clean), "\n\n")

cat("Distribución de CHOICE:\n")
print(table(sm_clean$choice))
cat("\n")

cat("Proporciones de CHOICE:\n")
print(round(prop.table(table(sm_clean$choice)), 4))
cat("\n")

# -----------------------------
# 5) Revisar rangos básicos
# -----------------------------
cat("=====================================================\n")
cat("Rangos de tiempos y costos\n")
cat("=====================================================\n")

rangos <- sm_clean |>
  summarise(
    train_tt_min = min(train_tt, na.rm = TRUE),
    train_tt_max = max(train_tt, na.rm = TRUE),
    train_co_min = min(train_co, na.rm = TRUE),
    train_co_max = max(train_co, na.rm = TRUE),
    
    sm_tt_min = min(sm_tt, na.rm = TRUE),
    sm_tt_max = max(sm_tt, na.rm = TRUE),
    sm_co_min = min(sm_co, na.rm = TRUE),
    sm_co_max = max(sm_co, na.rm = TRUE),
    
    car_tt_min = min(car_tt, na.rm = TRUE),
    car_tt_max = max(car_tt, na.rm = TRUE),
    car_co_min = min(car_co, na.rm = TRUE),
    car_co_max = max(car_co, na.rm = TRUE)
  )

print(rangos)
cat("\n")

# -----------------------------
# 6) Recodificaciones simples
# -----------------------------
# luggage en el archivo original suele venir como:
# 0 = none, 1 = one piece, 3 = several pieces
# Lo recodificamos a 0, 1, 2 para dejarlo más limpio
sm_clean <- sm_clean |>
  mutate(
    luggage = case_when(
      luggage == 0 ~ 0L,
      luggage == 1 ~ 1L,
      luggage == 3 ~ 2L,
      TRUE ~ NA_integer_
    )
  )

# -----------------------------
# 7) Selección de variables para modelamiento
# -----------------------------
# Dejamos una base limpia y reducida, centrada en lo que usarás en la tesis
sm_model <- sm_clean |>
  select(
    id, choice,
    train_av, car_av, sm_av,
    train_tt, train_co, train_he,
    sm_tt, sm_co, sm_he, sm_seats,
    car_tt, car_co,
    age, ga, luggage,
    male, income, purpose, first, group
  )

# -----------------------------
# 8) Chequeos finales
# -----------------------------
cat("=====================================================\n")
cat("Chequeos finales\n")
cat("=====================================================\n")
cat("Observaciones finales:", nrow(sm_model), "\n")
cat("Columnas finales:", ncol(sm_model), "\n\n")

cat("NA por variable en la base final:\n")
print(colSums(is.na(sm_model)))
cat("\n")

cat("Resumen de variables principales:\n")
print(summary(sm_model))
cat("\n")

# Aviso comparativo con el tamaño esperado según tesis/paper
if (nrow(sm_model) == 9036) {
  cat("OK: la muestra final coincide con 9036 observaciones.\n")
} else {
  cat("AVISO: la muestra final NO es 9036.\n")
  cat("Revisa filtros, nombres de variables o codificación del archivo.\n")
}
cat("\n")

# -----------------------------
# 9) Guardar resultados
# -----------------------------
write_csv(sm_model, archivo_salida_csv)
saveRDS(sm_model, archivo_salida_rds)

cat("=====================================================\n")
cat("Archivos guardados\n")
cat("=====================================================\n")
cat("CSV :", archivo_salida_csv, "\n")
cat("RDS :", archivo_salida_rds, "\n")
cat("Proceso terminado.\n")