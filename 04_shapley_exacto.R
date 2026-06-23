# =========================================================
# 04_shapley_exacto.R
# Shapley exacto global y local para el MNL base Swissmetro
# =========================================================

# 1. Cargar paquetes
library(apollo)
library(ggplot2)

set.seed(2026)

# 2. Definir rutas y opciones
calcular_local <- TRUE
hacer_graficos <- TRUE
reconstruir_cache <- FALSE
limpiar_archivos_apollo <- TRUE

archivo_entrada <- file.path("data", "processed", "swissmetro_clean.csv")
directorio_salida <- "output"
directorio_tablas <- file.path(directorio_salida, "tablas")
directorio_figuras <- file.path(directorio_salida, "figuras")
directorio_apollo <- file.path(directorio_salida, "apollo_coalitions")
directorio_cache <- file.path(directorio_salida, "cache_coaliciones")

dir.create(directorio_tablas, recursive = TRUE, showWarnings = FALSE)
dir.create(directorio_figuras, recursive = TRUE, showWarnings = FALSE)
dir.create(directorio_apollo, recursive = TRUE, showWarnings = FALSE)

if (reconstruir_cache && dir.exists(directorio_cache)) {
  unlink(directorio_cache, recursive = TRUE, force = TRUE)
}

dir.create(directorio_cache, recursive = TRUE, showWarnings = FALSE)

# 3. Definir jugadores y parametros del MNL
jugadores <- c("mode", "time", "cost", "age", "ga", "luggage", "seats")

mapa_parametros <- list(
  mode = c("asc_car", "asc_sm"),
  time = "b_time",
  cost = "b_cost",
  age = "b_age",
  ga = "b_ga",
  luggage = "b_luggage",
  seats = "b_seats"
)

todos_parametros <- unique(unlist(mapa_parametros, use.names = FALSE))

# Coeficientes del modelo completo usados como punto inicial.
# Estos valores reproducen la Tabla 7 de Salas et al. (2025).
coeficientes_iniciales <- c(
  asc_car = 2.132789,
  asc_sm = 1.870174,
  b_time = -0.012359,
  b_cost = -0.001011,
  b_age = 0.251225,
  b_ga = 6.656641,
  b_luggage = -0.083423,
  b_seats = 0.550311
)

# 4. Cargar datos depurados
df <- read.csv(archivo_entrada, stringsAsFactors = FALSE)
names(df) <- tolower(names(df))

df$id <- as.integer(df$id)
df$choice <- as.integer(df$choice)
df$train_av <- as.integer(df$train_av)
df$sm_av <- as.integer(df$sm_av)
df$car_av <- as.integer(df$car_av)
df$train_tt <- as.numeric(df$train_tt)
df$train_co <- as.numeric(df$train_co)
df$sm_tt <- as.numeric(df$sm_tt)
df$sm_co <- as.numeric(df$sm_co)
df$car_tt <- as.numeric(df$car_tt)
df$car_co <- as.numeric(df$car_co)
df$age <- as.numeric(df$age)
df$ga <- as.numeric(df$ga)
df$luggage <- as.numeric(df$luggage)
df$sm_seats <- as.numeric(df$sm_seats)

n_obs <- nrow(df)
ll_nulo <- n_obs * log(1 / 3)

# 5. Definir funciones auxiliares
completar_beta <- function(parametros_activos, valores_activos) {
  beta <- setNames(rep(0, length(todos_parametros)), todos_parametros)

  if (length(parametros_activos) > 0) {
    beta[parametros_activos] <- valores_activos
  }

  return(beta)
}

calcular_utilidades <- function(datos, beta) {
  V_train <-
    beta["b_time"] * datos$train_tt +
    beta["b_cost"] * datos$train_co +
    beta["b_age"] * datos$age +
    beta["b_ga"] * datos$ga

  V_sm <-
    beta["asc_sm"] +
    beta["b_time"] * datos$sm_tt +
    beta["b_cost"] * datos$sm_co +
    beta["b_ga"] * datos$ga +
    beta["b_seats"] * datos$sm_seats

  V_car <-
    beta["asc_car"] +
    beta["b_time"] * datos$car_tt +
    beta["b_cost"] * datos$car_co +
    beta["b_luggage"] * datos$luggage

  return(cbind(train = V_train, sm = V_sm, car = V_car))
}

calcular_probabilidades <- function(V) {
  V_centrado <- V - apply(V, 1, max)
  exp_V <- exp(V_centrado)
  return(exp_V / rowSums(exp_V))
}

calcular_loglik <- function(choice, probabilidades) {
  prob_elegida <- probabilidades[cbind(seq_along(choice), choice)]
  return(sum(log(pmax(prob_elegida, 1e-300))))
}

clave_coalicion <- function(S) {
  if (length(S) == 0) {
    return("EMPTY")
  }

  return(gsub("[^A-Za-z0-9_\\-]", "_", paste(sort(S), collapse = "__")))
}

parametros_coalicion <- function(S) {
  return(unique(unlist(mapa_parametros[S], use.names = FALSE)))
}

generar_coaliciones <- function(jugadores) {
  coaliciones <- list(character(0))

  for (k in seq_along(jugadores)) {
    coaliciones <- c(coaliciones, combn(jugadores, k, simplify = FALSE))
  }

  return(coaliciones)
}

limpiar_apollo <- function(nombre_modelo) {
  archivos <- list.files(
    directorio_apollo,
    pattern = paste0("^", nombre_modelo),
    full.names = TRUE
  )

  if (length(archivos) > 0) {
    unlink(archivos, force = TRUE)
  }
}

normalizar_resultado <- function(resultado) {
  if (is.null(resultado$coalicion) && !is.null(resultado$coalition)) {
    resultado$coalicion <- resultado$coalition
  }

  if (is.null(resultado$parametros_activos) && !is.null(resultado$active_params)) {
    resultado$parametros_activos <- resultado$active_params
  }

  if (is.null(resultado$probabilidades) && !is.null(resultado$probs)) {
    resultado$probabilidades <- resultado$probs
  }

  if (is.null(resultado$convergencia) && !is.null(resultado$converged)) {
    resultado$convergencia <- resultado$converged
  }

  if (is.null(resultado$metodo) && !is.null(resultado$method)) {
    resultado$metodo <- resultado$method
  }

  return(resultado)
}

# 6. Estimar modelos por coalicion
estimar_coalicion_vacia <- function(S, datos, archivo_cache) {
  beta <- completar_beta(character(0), numeric(0))
  probabilidades <- calcular_probabilidades(calcular_utilidades(datos, beta))
  ll_final <- calcular_loglik(datos$choice, probabilidades)

  resultado <- list(
    coalicion = S,
    key = clave_coalicion(S),
    parametros_activos = character(0),
    beta = beta,
    ll_final = ll_final,
    rho2 = 1 - ll_final / ll_nulo,
    probabilidades = probabilidades,
    convergencia = TRUE,
    metodo = "manual_empty"
  )

  saveRDS(resultado, archivo_cache)
  return(resultado)
}

estimar_coalicion <- function(S, datos) {
  key <- clave_coalicion(S)
  archivo_cache <- file.path(directorio_cache, paste0("coal_", key, ".rds"))

  if (file.exists(archivo_cache)) {
    return(normalizar_resultado(readRDS(archivo_cache)))
  }

  parametros_activos <- parametros_coalicion(S)

  if (length(parametros_activos) == 0) {
    return(estimar_coalicion_vacia(S, datos, archivo_cache))
  }

  database <- datos
  nombre_modelo <- paste0("smcoal_", key)

  apollo_control <- list(
    modelName = nombre_modelo,
    modelDescr = paste("Swissmetro coalition", key),
    indivID = "id",
    panelData = TRUE,
    outputDirectory = directorio_apollo
  )

  apollo_beta <- coeficientes_iniciales[parametros_activos]
  apollo_fixed <- c()

  apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))

    V <- list(train = 0, sm = 0, car = 0)

    if ("asc_sm" %in% names(apollo_beta)) {
      V[["sm"]] <- V[["sm"]] + asc_sm
    }

    if ("asc_car" %in% names(apollo_beta)) {
      V[["car"]] <- V[["car"]] + asc_car
    }

    if ("b_time" %in% names(apollo_beta)) {
      V[["train"]] <- V[["train"]] + b_time * train_tt
      V[["sm"]] <- V[["sm"]] + b_time * sm_tt
      V[["car"]] <- V[["car"]] + b_time * car_tt
    }

    if ("b_cost" %in% names(apollo_beta)) {
      V[["train"]] <- V[["train"]] + b_cost * train_co
      V[["sm"]] <- V[["sm"]] + b_cost * sm_co
      V[["car"]] <- V[["car"]] + b_cost * car_co
    }

    if ("b_age" %in% names(apollo_beta)) {
      V[["train"]] <- V[["train"]] + b_age * age
    }

    if ("b_ga" %in% names(apollo_beta)) {
      V[["train"]] <- V[["train"]] + b_ga * ga
      V[["sm"]] <- V[["sm"]] + b_ga * ga
    }

    if ("b_luggage" %in% names(apollo_beta)) {
      V[["car"]] <- V[["car"]] + b_luggage * luggage
    }

    if ("b_seats" %in% names(apollo_beta)) {
      V[["sm"]] <- V[["sm"]] + b_seats * sm_seats
    }

    mnl_settings <- list(
      alternatives = c(train = 1, sm = 2, car = 3),
      avail = list(train = train_av, sm = sm_av, car = car_av),
      choiceVar = choice,
      V = V
    )

    P <- list()
    P[["model"]] <- apollo_mnl(mnl_settings, functionality)
    P <- apollo_panelProd(P, apollo_inputs, functionality)
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
  }

  invisible(capture.output({
    apollo_inputs <- apollo_validateInputs()
  }))

  invisible(capture.output({
    modelo <- apollo_estimate(
      apollo_beta,
      apollo_fixed,
      apollo_probabilities,
      apollo_inputs
    )
  }))

  beta <- completar_beta(parametros_activos, modelo$estimate)
  probabilidades <- calcular_probabilidades(calcular_utilidades(datos, beta))
  ll_final <- calcular_loglik(datos$choice, probabilidades)

  resultado <- list(
    coalicion = S,
    key = key,
    parametros_activos = parametros_activos,
    beta = beta,
    ll_final = ll_final,
    rho2 = 1 - ll_final / ll_nulo,
    probabilidades = probabilidades,
    convergencia = TRUE,
    metodo = "Apollo_warmstart"
  )

  saveRDS(resultado, archivo_cache)

  if (limpiar_archivos_apollo) {
    limpiar_apollo(nombre_modelo)
  }

  return(resultado)
}

# 7. Calcular Shapley exacto
peso_shapley <- function(s, p) {
  return(factorial(s) * factorial(p - s - 1) / factorial(p))
}

calcular_shapley_global <- function(resultados, coaliciones, jugadores) {
  p <- length(jugadores)
  phi <- setNames(rep(0, p), jugadores)

  for (jugador in jugadores) {
    coaliciones_sin_jugador <- coaliciones[
      !sapply(coaliciones, function(S) jugador %in% S)
    ]

    for (S in coaliciones_sin_jugador) {
      key_S <- clave_coalicion(S)
      key_S_mas_jugador <- clave_coalicion(c(S, jugador))
      w <- peso_shapley(length(S), p)

      phi[jugador] <- phi[jugador] +
        w * (resultados[[key_S_mas_jugador]]$rho2 - resultados[[key_S]]$rho2)
    }
  }

  tabla <- data.frame(
    player = names(phi),
    shapley_value = as.numeric(phi),
    share_percent = 100 * as.numeric(phi) / sum(phi),
    row.names = NULL
  )

  tabla <- tabla[order(-tabla$shapley_value), ]
  return(tabla)
}

calcular_shapley_local <- function(resultados, coaliciones, jugadores, n_obs) {
  p <- length(jugadores)
  phi <- array(
    0,
    dim = c(n_obs, 3, p),
    dimnames = list(NULL, c("train", "sm", "car"), jugadores)
  )

  for (i in seq_along(jugadores)) {
    jugador <- jugadores[i]
    cat("Calculando SHAP local para:", jugador, "\n")

    coaliciones_sin_jugador <- coaliciones[
      !sapply(coaliciones, function(S) jugador %in% S)
    ]

    for (S in coaliciones_sin_jugador) {
      key_S <- clave_coalicion(S)
      key_S_mas_jugador <- clave_coalicion(c(S, jugador))
      w <- peso_shapley(length(S), p)

      phi[, , i] <- phi[, , i] +
        w * (
          resultados[[key_S_mas_jugador]]$probabilidades -
            resultados[[key_S]]$probabilidades
        )
    }
  }

  return(phi)
}

# 8. Preparar tablas, chequeos y graficos
resumir_coaliciones <- function(resultados) {
  tabla <- do.call(
    rbind,
    lapply(resultados, function(x) {
      data.frame(
        key = x$key,
        coalition_size = length(x$coalicion),
        coalition = if (length(x$coalicion) == 0) "EMPTY" else paste(x$coalicion, collapse = ", "),
        ll_final = x$ll_final,
        rho2 = x$rho2,
        converged = x$convergencia,
        method = x$metodo,
        stringsAsFactors = FALSE
      )
    })
  )

  tabla <- tabla[order(tabla$coalition_size, tabla$key), ]
  row.names(tabla) <- NULL
  return(tabla)
}

comparar_modelo_completo <- function(resultados, datos) {
  key_completo <- clave_coalicion(jugadores)
  modelo_completo <- resultados[[key_completo]]

  probabilidades_benchmark <- calcular_probabilidades(
    calcular_utilidades(datos, coeficientes_iniciales)
  )
  ll_benchmark <- calcular_loglik(datos$choice, probabilidades_benchmark)
  rho2_benchmark <- 1 - ll_benchmark / ll_nulo

  parametros <- data.frame(
    parametro = names(coeficientes_iniciales),
    benchmark = as.numeric(coeficientes_iniciales),
    estimate = as.numeric(modelo_completo$beta[names(coeficientes_iniciales)]),
    stringsAsFactors = FALSE
  )
  parametros$diff <- parametros$estimate - parametros$benchmark

  metricas <- data.frame(
    metric = c("ll_benchmark", "ll_full", "rho2_benchmark", "rho2_full"),
    value = c(ll_benchmark, modelo_completo$ll_final, rho2_benchmark, modelo_completo$rho2),
    stringsAsFactors = FALSE
  )

  return(list(parametros = parametros, metricas = metricas, modelo = modelo_completo))
}

chequear_eficiencia <- function(resultados, shapley_global, shapley_local) {
  key_completo <- clave_coalicion(jugadores)

  error_global <- sum(shapley_global$shapley_value) -
    (resultados[[key_completo]]$rho2 - resultados[["EMPTY"]]$rho2)

  chequeos <- data.frame(
    check = "global_efficiency_error",
    value = error_global,
    stringsAsFactors = FALSE
  )

  if (!is.null(shapley_local)) {
    suma_local <- apply(shapley_local, c(1, 2), sum)
    objetivo_local <- resultados[[key_completo]]$probabilidades -
      resultados[["EMPTY"]]$probabilidades

    chequeos <- rbind(
      chequeos,
      data.frame(
        check = "max_local_efficiency_error",
        value = max(abs(suma_local - objetivo_local)),
        stringsAsFactors = FALSE
      )
    )
  }

  return(chequeos)
}

preparar_datos_locales <- function(datos, shapley_local) {
  tabla <- rbind(
    data.frame(alternative = "car", feature = "CAR_CO", value = datos$car_co, shap = shapley_local[, "car", "cost"]),
    data.frame(alternative = "car", feature = "CAR_TT", value = datos$car_tt, shap = shapley_local[, "car", "time"]),
    data.frame(alternative = "car", feature = "LUGGAGE", value = datos$luggage, shap = shapley_local[, "car", "luggage"]),
    data.frame(alternative = "train", feature = "GA", value = datos$ga, shap = shapley_local[, "train", "ga"]),
    data.frame(alternative = "train", feature = "TRAIN_TT", value = datos$train_tt, shap = shapley_local[, "train", "time"]),
    data.frame(alternative = "train", feature = "TRAIN_CO", value = datos$train_co, shap = shapley_local[, "train", "cost"]),
    data.frame(alternative = "train", feature = "AGE", value = datos$age, shap = shapley_local[, "train", "age"]),
    data.frame(alternative = "sm", feature = "SM_TT", value = datos$sm_tt, shap = shapley_local[, "sm", "time"]),
    data.frame(alternative = "sm", feature = "GA", value = datos$ga, shap = shapley_local[, "sm", "ga"]),
    data.frame(alternative = "sm", feature = "SM_CO", value = datos$sm_co, shap = shapley_local[, "sm", "cost"]),
    data.frame(alternative = "sm", feature = "SM_SEATS", value = datos$sm_seats, shap = shapley_local[, "sm", "seats"])
  )

  tabla$mean_abs_shap <- ave(abs(tabla$shap), tabla$alternative, tabla$feature, FUN = mean)
  return(tabla)
}

guardar_graficos <- function(shapley_global, datos_locales) {
  grafico_global <- ggplot(
    shapley_global,
    aes(x = reorder(player, shapley_value), y = share_percent)
  ) +
    geom_col(fill = "#3366A8", width = 0.7) +
    coord_flip() +
    labs(
      title = "Shapley global exacto",
      x = NULL,
      y = "Participacion (%)"
    ) +
    theme_minimal()

  ggsave(
    file.path(directorio_figuras, "swissmetro_shapley_global_exact.png"),
    grafico_global,
    width = 7,
    height = 4.5,
    dpi = 300
  )

  if (!is.null(datos_locales)) {
    datos_locales$feature <- reorder(datos_locales$feature, datos_locales$mean_abs_shap)

    grafico_local <- ggplot(
      datos_locales,
      aes(x = shap, y = feature, color = value)
    ) +
      geom_point(
        alpha = 0.45,
        size = 0.6,
        position = position_jitter(height = 0.18, width = 0)
      ) +
      facet_wrap(~alternative, scales = "free_y") +
      scale_color_gradient(low = "#3366A8", high = "#B33A3A") +
      labs(
        title = "SHAP local exacto por alternativa",
        x = "Contribucion a la probabilidad predicha",
        y = NULL,
        color = "Valor"
      ) +
      theme_minimal()

    ggsave(
      file.path(directorio_figuras, "swissmetro_shapley_local_summary.png"),
      grafico_local,
      width = 10,
      height = 6,
      dpi = 300
    )
  }
}

# 9. Estimar todas las coaliciones
apollo_initialise()

coaliciones <- generar_coaliciones(jugadores)
resultados <- vector("list", length(coaliciones))
names(resultados) <- sapply(coaliciones, clave_coalicion)

cat("=====================================================\n")
cat("Swissmetro: Shapley exacto global y local\n")
cat("=====================================================\n")
cat("Observaciones:", n_obs, "\n")
cat("IDs unicos:", length(unique(df$id)), "\n")
cat("Jugadores:", paste(jugadores, collapse = ", "), "\n")
cat("Coaliciones:", length(coaliciones), "de", 2^length(jugadores), "\n\n")

for (i in seq_along(coaliciones)) {
  key <- names(resultados)[i]
  cat(sprintf("[%3d/%3d] Coalicion: %s\n", i, length(coaliciones), key))
  resultados[[key]] <- estimar_coalicion(coaliciones[[i]], df)
}

# 10. Calcular resultados exactos
resumen_coaliciones <- resumir_coaliciones(resultados)
comparacion_completa <- comparar_modelo_completo(resultados, df)
shapley_global <- calcular_shapley_global(resultados, coaliciones, jugadores)

shapley_local <- NULL
datos_locales <- NULL

if (calcular_local) {
  shapley_local <- calcular_shapley_local(resultados, coaliciones, jugadores, n_obs)
  datos_locales <- preparar_datos_locales(df, shapley_local)
}

chequeos <- chequear_eficiencia(resultados, shapley_global, shapley_local)

# La eficiencia es la propiedad central del calculo exacto.
if (max(abs(chequeos$value)) > 1e-8) {
  stop("Fallo el chequeo de eficiencia Shapley.", call. = FALSE)
}

# 11. Guardar tablas y objetos
write.csv(
  resumen_coaliciones,
  file.path(directorio_tablas, "swissmetro_coalitions_summary.csv"),
  row.names = FALSE
)

write.csv(
  comparacion_completa$parametros,
  file.path(directorio_tablas, "swissmetro_full_model_vs_benchmark.csv"),
  row.names = FALSE
)

write.csv(
  comparacion_completa$metricas,
  file.path(directorio_tablas, "swissmetro_full_model_metrics.csv"),
  row.names = FALSE
)

write.csv(
  shapley_global,
  file.path(directorio_tablas, "swissmetro_shapley_global_exact.csv"),
  row.names = FALSE
)

write.csv(
  chequeos,
  file.path(directorio_tablas, "swissmetro_shapley_exact_checks.csv"),
  row.names = FALSE
)

if (!is.null(datos_locales)) {
  write.csv(
    datos_locales,
    file.path(directorio_tablas, "swissmetro_shapley_local_plot_data.csv"),
    row.names = FALSE
  )

  saveRDS(
    shapley_local,
    file.path(directorio_salida, "swissmetro_shapley_local_exact.rds")
  )
}

objeto_final <- list(
  players = jugadores,
  param_map = mapa_parametros,
  benchmark_start = coeficientes_iniciales,
  ll_null = ll_nulo,
  coalition_summary = resumen_coaliciones,
  full_model = comparacion_completa$modelo,
  full_model_comparison = comparacion_completa$parametros,
  full_model_metrics = comparacion_completa$metricas,
  shapley_global = shapley_global,
  shapley_local = shapley_local,
  checks = chequeos
)

saveRDS(
  objeto_final,
  file.path(directorio_salida, "swissmetro_shapley_exact_results.rds")
)

# 12. Guardar graficos
if (hacer_graficos) {
  guardar_graficos(shapley_global, datos_locales)
}

# 13. Imprimir resultados principales
cat("\n=====================================================\n")
cat("Shapley exacto terminado\n")
cat("=====================================================\n")
cat("Observaciones usadas:", n_obs, "\n")
cat("Error eficiencia global:", chequeos$value[chequeos$check == "global_efficiency_error"], "\n")

if ("max_local_efficiency_error" %in% chequeos$check) {
  cat("Error maximo eficiencia local:", chequeos$value[chequeos$check == "max_local_efficiency_error"], "\n")
}

cat("\nComparacion modelo completo con Salas et al. (2025):\n")
print(comparacion_completa$parametros)
cat("\n")

cat("Shapley global exacto:\n")
print(shapley_global)
cat("\n")

cat("Tablas guardadas en:", directorio_tablas, "\n")
cat("Graficos guardados en:", directorio_figuras, "\n")
