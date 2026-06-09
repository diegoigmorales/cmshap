# =========================================================
# 03_reproducir_salas_swissmetro.R
# Reproducción del MNL base de Salas et al. (2025) en Swissmetro
# =========================================================

# -----------------------------
# 0) Paquetes
# -----------------------------
library(apollo)
library(dplyr)
library(readr)

# -----------------------------
# 1) Cargar base limpia
# -----------------------------
sm <- readRDS("data/processed/swissmetro_clean.rds")

cat("=====================================================\n")
cat("Swissmetro: reproducción MNL base\n")
cat("=====================================================\n")
cat("Observaciones:", nrow(sm), "\n")
cat("IDs únicos:", dplyr::n_distinct(sm$id), "\n\n")

# Verificación rápida
cat("Distribución de CHOICE:\n")
print(table(sm$choice))
cat("\n")

# -----------------------------
# 2) Preparar datos
# -----------------------------
# Nos aseguramos de que las variables estén como numéricas
sm <- sm |>
  mutate(
    choice    = as.integer(choice),
    id        = as.integer(id),
    train_av  = as.integer(train_av),
    sm_av     = as.integer(sm_av),
    car_av    = as.integer(car_av),
    train_tt  = as.numeric(train_tt),
    train_co  = as.numeric(train_co),
    sm_tt     = as.numeric(sm_tt),
    sm_co     = as.numeric(sm_co),
    car_tt    = as.numeric(car_tt),
    car_co    = as.numeric(car_co),
    age       = as.numeric(age),
    ga        = as.numeric(ga),
    luggage   = as.numeric(luggage),
    sm_seats  = as.numeric(sm_seats)
  )

database <- sm

# -----------------------------
# 3) Inicializar Apollo
# -----------------------------
apollo_initialise()

apollo_control <- list(
  modelName       = "Swissmetro_MNL_Salas2025",
  modelDescr      = "Reproducción MNL base de Salas et al. (2025) en Swissmetro",
  indivID         = "id",
  outputDirectory = "output"
)

# -----------------------------
# 4) Parámetros
# -----------------------------
# Train queda como alternativa base: no lleva ASC
apollo_beta <- c(
  asc_car     = 0,
  asc_sm      = 0,
  b_time      = 0,
  b_cost      = 0,
  b_age       = 0,
  b_ga        = 0,
  b_luggage   = 0,
  b_seats     = 0
)

apollo_fixed <- c()

# -----------------------------
# 5) Validar inputs
# -----------------------------
apollo_inputs <- apollo_validateInputs()

# -----------------------------
# 6) Definir probabilidades
# -----------------------------
apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  P <- list()
  
  # Utilidades según Table 5 de Salas et al. (2025)
  # Car, Train, Swissmetro
  V <- list()
  
  V[["train"]] <- 
    b_time * train_tt +
    b_cost * train_co +
    b_age  * age +
    b_ga   * ga
  
  V[["sm"]] <- 
    asc_sm +
    b_time  * sm_tt +
    b_cost  * sm_co +
    b_ga    * ga +
    b_seats * sm_seats
  
  V[["car"]] <- 
    asc_car +
    b_time    * car_tt +
    b_cost    * car_co +
    b_luggage * luggage
  
  mnl_settings <- list(
    alternatives = c(train = 1, sm = 2, car = 3),
    avail        = list(train = train_av, sm = sm_av, car = car_av),
    choiceVar    = choice,
    V            = V
  )
  
  P[["model"]] <- apollo_mnl(mnl_settings, functionality)
  
  # Swissmetro tiene varias observaciones por individuo
  P <- apollo_panelProd(P, apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  
  return(P)
}

# -----------------------------
# 7) Estimar modelo
# -----------------------------
modelo <- apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs
)

# -----------------------------
# 8) Mostrar y guardar salida
# -----------------------------
apollo_modelOutput(modelo)
apollo_saveOutput(modelo)

# -----------------------------
# 9) Extraer resultados principales
# -----------------------------
coef_est <- modelo$estimate
se_est   <- modelo$se
t_est    <- coef_est / se_est

resultados <- data.frame(
  parametro = names(coef_est),
  estimate  = as.numeric(coef_est),
  se        = as.numeric(se_est),
  t_value   = as.numeric(t_est)
)

cat("=====================================================\n")
cat("Coeficientes estimados\n")
cat("=====================================================\n")
print(resultados)
cat("\n")

# -----------------------------
# 10) Log-likelihood y pseudo-R2
# -----------------------------
# Como después de depurar todas las observaciones tienen las 3 alternativas disponibles,
# el log-likelihood "at zero" con todas las utilidades en 0 es:
# LL0 = N * log(1/3)
ll_zero  <- nrow(sm) * log(1/3)

# En Apollo, el máximo suele quedar en:
ll_final <- modelo$maximum

pseudo_r2 <- 1 - (ll_final / ll_zero)

cat("=====================================================\n")
cat("Medidas de ajuste\n")
cat("=====================================================\n")
cat("LL(0)          =", round(ll_zero, 3), "\n")
cat("LL(final)      =", round(ll_final, 3), "\n")
cat("Pseudo-R2 MF   =", round(pseudo_r2, 3), "\n\n")

# -----------------------------
# 11) Comparación con Table 7 de Salas et al. (2025)
# -----------------------------
objetivo <- data.frame(
  parametro = c("asc_car", "asc_sm", "b_cost", "b_time", "b_age", "b_luggage", "b_ga", "b_seats"),
  target    = c(2.137, 1.871, -0.001, -0.012, 0.252, -0.087, 6.650, 0.547)
)

comparacion <- resultados |>
  select(parametro, estimate) |>
  left_join(objetivo, by = "parametro") |>
  mutate(
    diferencia = estimate - target
  )

cat("=====================================================\n")
cat("Comparación con Table 7\n")
cat("=====================================================\n")
print(comparacion)
cat("\n")

# -----------------------------
# 12) Guardar tablas
# -----------------------------
dir.create("output/tablas", recursive = TRUE, showWarnings = FALSE)

write_csv(resultados,   "output/tablas/swissmetro_mnl_coeficientes.csv")
write_csv(comparacion,  "output/tablas/swissmetro_mnl_comparacion_table7.csv")

cat("Archivos guardados en output/tablas/\n")
cat("Fin de la reproducción del MNL base.\n")