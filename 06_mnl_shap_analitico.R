# =========================================================
# 06_mnl_shap_analitico.R
# MNL-SHAP analitico para Swissmetro sin coaliciones
# =========================================================

# 1. Cargar paquetes
library(apollo)
library(ggplot2)

set.seed(2026)

# 2. Definir rutas
archivo_entrada <- "swissmetro_limpio.csv"
directorio_salida <- "output"
directorio_tablas <- file.path(directorio_salida, "tablas")
directorio_figuras <- file.path(directorio_salida, "figuras")
directorio_apollo <- file.path(directorio_salida, "apollo_mnl_shap_analitico")

dir.create(directorio_tablas, recursive = TRUE, showWarnings = FALSE)
dir.create(directorio_figuras, recursive = TRUE, showWarnings = FALSE)
dir.create(directorio_apollo, recursive = TRUE, showWarnings = FALSE)

# 3. Cargar datos resultantes de 02_depurar.R
if (!file.exists(archivo_entrada)) {
  source("02_depurar.R")
}

df <- read.csv(archivo_entrada, stringsAsFactors = FALSE)

df$ID <- seq_len(nrow(df))
df$TRAIN_AV <- 1
df$SM_AV <- 1
df$CAR_AV <- 1

database <- df

# 4. Especificar modelo MNL en Apollo
apollo_initialise()

apollo_control <- list(
  modelName = "Swissmetro_MNL_SHAP_Analitico",
  modelDescr = "MNL base Swissmetro para MNL-SHAP analitico",
  indivID = "ID",
  outputDirectory = directorio_apollo
)

apollo_beta <- c(
  asc_car = 0,
  asc_sm = 0,
  b_time = 0,
  b_cost = 0,
  b_age = 0,
  b_ga = 0,
  b_luggage = 0,
  b_seats = 0
)

apollo_fixed <- c()

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  V <- list()

  V[["train"]] <-
    b_time * TRAIN_TT +
    b_cost * TRAIN_CO +
    b_age * AGE +
    b_ga * GA

  V[["sm"]] <-
    asc_sm +
    b_time * SM_TT +
    b_cost * SM_CO +
    b_ga * GA +
    b_seats * SM_SEATS

  V[["car"]] <-
    asc_car +
    b_time * CAR_TT +
    b_cost * CAR_CO +
    b_luggage * LUGGAGE

  mnl_settings <- list(
    alternatives = c(train = 1, sm = 2, car = 3),
    avail = list(train = TRAIN_AV, sm = SM_AV, car = CAR_AV),
    choiceVar = CHOICE,
    V = V
  )

  P <- list()
  P[["model"]] <- apollo_mnl(mnl_settings, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

apollo_inputs <- apollo_validateInputs()

# 5. Estimar modelo MNL
modelo <- apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs
)

apollo_modelOutput(modelo)
apollo_saveOutput(modelo)

coeficientes <- modelo$estimate

tabla_coeficientes <- data.frame(
  parametro = names(coeficientes),
  estimate = as.numeric(coeficientes),
  row.names = NULL
)

if (!is.null(modelo$se)) {
  tabla_coeficientes$se <- as.numeric(modelo$se[names(coeficientes)])
  tabla_coeficientes$t_value <- tabla_coeficientes$estimate / tabla_coeficientes$se
}

# 6. Construir utilidades
alternativas <- c("train", "sm", "car")

V <- cbind(
  train =
    coeficientes["b_time"] * df$TRAIN_TT +
    coeficientes["b_cost"] * df$TRAIN_CO +
    coeficientes["b_age"] * df$AGE +
    coeficientes["b_ga"] * df$GA,
  sm =
    coeficientes["asc_sm"] +
    coeficientes["b_time"] * df$SM_TT +
    coeficientes["b_cost"] * df$SM_CO +
    coeficientes["b_ga"] * df$GA +
    coeficientes["b_seats"] * df$SM_SEATS,
  car =
    coeficientes["asc_car"] +
    coeficientes["b_time"] * df$CAR_TT +
    coeficientes["b_cost"] * df$CAR_CO +
    coeficientes["b_luggage"] * df$LUGGAGE
)

variables_shap <- data.frame(
  atributo = c(
    "TRAIN_TT", "SM_TT", "CAR_TT",
    "TRAIN_CO", "SM_CO", "CAR_CO",
    "AGE_train", "GA_train", "GA_sm",
    "LUGGAGE_car", "SEATS_sm"
  ),
  grupo = c(
    "TT", "TT", "TT",
    "COST", "COST", "COST",
    "AGE", "GA", "GA",
    "LUGGAGE", "SEATS"
  ),
  alternativa = c(
    "train", "sm", "car",
    "train", "sm", "car",
    "train", "train", "sm",
    "car", "sm"
  ),
  variable = c(
    "TRAIN_TT", "SM_TT", "CAR_TT",
    "TRAIN_CO", "SM_CO", "CAR_CO",
    "AGE", "GA", "GA",
    "LUGGAGE", "SM_SEATS"
  ),
  parametro = c(
    "b_time", "b_time", "b_time",
    "b_cost", "b_cost", "b_cost",
    "b_age", "b_ga", "b_ga",
    "b_luggage", "b_seats"
  ),
  stringsAsFactors = FALSE
)

n_obs <- nrow(df)
n_alt <- length(alternativas)
n_var <- nrow(variables_shap)

X <- sapply(variables_shap$variable, function(x) df[[x]])
colnames(X) <- variables_shap$atributo

media_X <- colMeans(X)
X_centrado <- sweep(X, 2, media_X, "-")

beta <- matrix(
  0,
  nrow = n_alt,
  ncol = n_var,
  dimnames = list(alternativas, variables_shap$atributo)
)

for (i in seq_len(n_var)) {
  beta[variables_shap$alternativa[i], variables_shap$atributo[i]] <-
    coeficientes[variables_shap$parametro[i]]
}

ASC <- c(
  train = 0,
  sm = coeficientes["asc_sm"],
  car = coeficientes["asc_car"]
)

V0 <- as.numeric(ASC + beta %*% media_X)
names(V0) <- alternativas

# 7. Calcular Utility-SHAP
phi_V <- array(
  0,
  dim = c(n_obs, n_alt, n_var),
  dimnames = list(NULL, alternativas, variables_shap$atributo)
)

for (j in alternativas) {
  phi_V[, j, ] <- sweep(X_centrado, 2, beta[j, ], "*")
}

V_reconstruido <- sweep(apply(phi_V, c(1, 2), sum), 2, V0, "+")
error_utilidad <- max(abs(V - V_reconstruido))

# 8. Calcular probabilidades MNL
V_centrado <- V - apply(V, 1, max)
exp_V <- exp(V_centrado)
P_mnl <- exp_V / rowSums(exp_V)

# 9. Calcular Probability-SHAP
phi_P <- array(
  0,
  dim = c(n_obs, n_alt, n_var),
  dimnames = list(NULL, alternativas, variables_shap$atributo)
)

for (m in seq_len(n_var)) {
  beta_promedio <- as.numeric(P_mnl %*% beta[, m])

  for (j in alternativas) {
    phi_P[, j, m] <-
      X_centrado[, m] *
      P_mnl[, j] *
      (beta[j, m] - beta_promedio)
  }
}

# 10. Calcular importancia global
calcular_importancia <- function(phi) {
  importancia <- data.frame(
    atributo = variables_shap$atributo,
    grupo = variables_shap$grupo,
    alternativa = variables_shap$alternativa,
    variable = variables_shap$variable,
    parametro = variables_shap$parametro,
    beta = NA_real_,
    media_base = media_X,
    importancia_global = NA_real_,
    importancia_alternativa = NA_real_,
    stringsAsFactors = FALSE
  )

  for (m in seq_len(n_var)) {
    importancia$beta[m] <- beta[variables_shap$alternativa[m], m]
    importancia$importancia_global[m] <- mean(abs(phi[, , m]))
    importancia$importancia_alternativa[m] <- mean(abs(phi[, variables_shap$alternativa[m], m]))
  }

  importancia <- importancia[order(-importancia$importancia_global), ]
  importancia$ranking <- seq_len(nrow(importancia))
  importancia <- importancia[
    c(
      "ranking", "atributo", "grupo", "alternativa", "variable",
      "parametro", "beta", "media_base",
      "importancia_global", "importancia_alternativa"
    )
  ]

  return(importancia)
}

importancia_utility <- calcular_importancia(phi_V)
importancia_probability <- calcular_importancia(phi_P)

diagnosticos <- data.frame(
  indicador = c(
    "observaciones",
    "variables_shap",
    "error_reconstruccion_utilidad",
    "log_likelihood_final"
  ),
  valor = c(
    n_obs,
    n_var,
    error_utilidad,
    modelo$maximum
  )
)

# 11. Exportar tablas
write.csv(
  tabla_coeficientes,
  file.path(directorio_tablas, "mnl_shap_analitico_coeficientes.csv"),
  row.names = FALSE
)

write.csv(
  diagnosticos,
  file.path(directorio_tablas, "mnl_shap_analitico_diagnosticos.csv"),
  row.names = FALSE
)

write.csv(
  importancia_utility,
  file.path(directorio_tablas, "mnl_shap_analitico_ranking_utility.csv"),
  row.names = FALSE
)

write.csv(
  importancia_probability,
  file.path(directorio_tablas, "mnl_shap_analitico_ranking_probability.csv"),
  row.names = FALSE
)

# 12. Preparar datos para graficos tipo jitter
datos_phi_V <- expand.grid(
  observacion = seq_len(n_obs),
  alternativa = alternativas,
  atributo = variables_shap$atributo,
  stringsAsFactors = FALSE
)

datos_phi_V$phi <- as.vector(phi_V)
datos_phi_V <- merge(datos_phi_V, variables_shap, by = "atributo")
datos_phi_V <- datos_phi_V[datos_phi_V$alternativa.x == datos_phi_V$alternativa.y, ]
names(datos_phi_V)[names(datos_phi_V) == "alternativa.x"] <- "alternativa"

datos_phi_P <- expand.grid(
  observacion = seq_len(n_obs),
  alternativa = alternativas,
  atributo = variables_shap$atributo,
  stringsAsFactors = FALSE
)

datos_phi_P$phi <- as.vector(phi_P)
datos_phi_P <- merge(datos_phi_P, variables_shap, by = "atributo")
names(datos_phi_P)[names(datos_phi_P) == "alternativa.x"] <- "alternativa"

if (nrow(datos_phi_V) > 60000) {
  datos_phi_V <- datos_phi_V[sample(seq_len(nrow(datos_phi_V)), 60000), ]
}

if (nrow(datos_phi_P) > 60000) {
  datos_phi_P <- datos_phi_P[sample(seq_len(nrow(datos_phi_P)), 60000), ]
}

write.csv(
  datos_phi_V,
  file.path(directorio_tablas, "mnl_shap_analitico_plotdata_utility.csv"),
  row.names = FALSE
)

write.csv(
  datos_phi_P,
  file.path(directorio_tablas, "mnl_shap_analitico_plotdata_probability.csv"),
  row.names = FALSE
)

# 13. Guardar graficos de importancia global
grafico_utility <- ggplot(
  importancia_utility,
  aes(x = reorder(atributo, importancia_global), y = importancia_global, fill = grupo)
) +
  geom_col(width = 0.75) +
  coord_flip() +
  labs(
    x = NULL,
    y = "mean(abs(phi^V))",
    fill = "Grupo",
    title = "Importancia global Utility-SHAP"
  ) +
  theme_minimal()

ggsave(
  file.path(directorio_figuras, "mnl_shap_analitico_importancia_utility.png"),
  grafico_utility,
  width = 8.5,
  height = 5.5,
  dpi = 300
)

grafico_probability <- ggplot(
  importancia_probability,
  aes(x = reorder(atributo, importancia_global), y = importancia_global, fill = grupo)
) +
  geom_col(width = 0.75) +
  coord_flip() +
  labs(
    x = NULL,
    y = "mean(abs(phi^P))",
    fill = "Grupo",
    title = "Importancia global Probability-SHAP"
  ) +
  theme_minimal()

ggsave(
  file.path(directorio_figuras, "mnl_shap_analitico_importancia_probability.png"),
  grafico_probability,
  width = 8.5,
  height = 5.5,
  dpi = 300
)

# 14. Guardar graficos tipo jitter
datos_phi_V$atributo <- factor(
  datos_phi_V$atributo,
  levels = rev(importancia_utility$atributo)
)

grafico_jitter_utility <- ggplot(
  datos_phi_V,
  aes(x = atributo, y = phi, color = alternativa)
) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.25) +
  geom_jitter(width = 0.22, height = 0, alpha = 0.24, size = 0.45) +
  coord_flip() +
  labs(
    x = NULL,
    y = "phi^V",
    color = "Alternativa",
    title = "Jitter Utility-SHAP"
  ) +
  theme_minimal()

ggsave(
  file.path(directorio_figuras, "mnl_shap_analitico_jitter_utility.png"),
  grafico_jitter_utility,
  width = 9,
  height = 6,
  dpi = 300
)

datos_phi_P$atributo <- factor(
  datos_phi_P$atributo,
  levels = rev(importancia_probability$atributo)
)

grafico_jitter_probability <- ggplot(
  datos_phi_P,
  aes(x = atributo, y = phi, color = alternativa)
) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.25) +
  geom_jitter(width = 0.22, height = 0, alpha = 0.24, size = 0.45) +
  coord_flip() +
  labs(
    x = NULL,
    y = "phi^P",
    color = "Alternativa",
    title = "Jitter Probability-SHAP"
  ) +
  theme_minimal()

ggsave(
  file.path(directorio_figuras, "mnl_shap_analitico_jitter_probability.png"),
  grafico_jitter_probability,
  width = 9,
  height = 6,
  dpi = 300
)

# 15. Imprimir resultados principales
cat("=====================================================\n")
cat("MNL-SHAP analitico terminado\n")
cat("=====================================================\n")
cat("Observaciones usadas:", n_obs, "\n")
cat("Error maximo de reconstruccion de utilidad:", error_utilidad, "\n\n")

cat("Ranking Utility-SHAP:\n")
print(importancia_utility)
cat("\n")

cat("Ranking Probability-SHAP:\n")
print(importancia_probability)
cat("\n")

cat("Tablas guardadas en:", directorio_tablas, "\n")
cat("Graficos guardados en:", directorio_figuras, "\n")
