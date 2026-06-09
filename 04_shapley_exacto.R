# =========================================================
# Swissmetro: Shapley exacto global + local con Apollo
# Warm starts + gráficos robustos (sin shapviz)
# Incluye una curva de predicción usando apollo_prediction()
# con fallback manual si Apollo no devuelve el formato esperado
# =========================================================

# -----------------------------
# 0) Paquetes
# -----------------------------
req_pkgs <- c(
  "apollo", "dplyr", "readr", "tidyr",
  "ggplot2", "ggbeeswarm", "patchwork", "scales"
)

miss_pkgs <- req_pkgs[!sapply(req_pkgs, requireNamespace, quietly = TRUE)]
if (length(miss_pkgs) > 0) {
  stop(
    paste0(
      "Faltan paquetes: ",
      paste(miss_pkgs, collapse = ", "),
      ". Instálalos antes de correr el script."
    )
  )
}

library(apollo)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(patchwork)
library(scales)

# -----------------------------
# 1) Cargar base limpia
# -----------------------------
sm <- readRDS("data/processed/swissmetro_clean.rds")

sm <- sm |>
  mutate(
    id        = as.integer(id),
    choice    = as.integer(choice),   # 1=train, 2=sm, 3=car
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

N <- nrow(sm)

cat("=====================================================\n")
cat("Swissmetro: exacto global + local con Apollo y warm starts\n")
cat("=====================================================\n")
cat("Observaciones:", N, "\n")
cat("IDs únicos:", dplyr::n_distinct(sm$id), "\n\n")

# -----------------------------
# 2) Configuración
# -----------------------------
players <- c("mode", "time", "cost", "age", "ga", "luggage", "seats")
P <- length(players)

ll_null <- N * log(1 / 3)

# Benchmark validado previamente
benchmark_start <- c(
  asc_car    =  2.132789,
  asc_sm     =  1.870174,
  b_time     = -0.012359,
  b_cost     = -0.001011,
  b_age      =  0.251225,
  b_ga       =  6.656641,
  b_luggage  = -0.083423,
  b_seats    =  0.550311
)

dir.create("output", showWarnings = FALSE)
dir.create("output/tablas", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figuras", recursive = TRUE, showWarnings = FALSE)
dir.create("output/apollo_coalitions", recursive = TRUE, showWarnings = FALSE)

unlink("output/cache_coaliciones", recursive = TRUE, force = TRUE)
dir.create("output/cache_coaliciones", recursive = TRUE, showWarnings = FALSE)

CALC_LOCAL <- TRUE
MAKE_PLOTS <- TRUE

apollo_initialise()

# -----------------------------
# 3) Funciones auxiliares
# -----------------------------
all_param_names <- c(
  "asc_car", "asc_sm",
  "b_time", "b_cost", "b_age", "b_ga", "b_luggage", "b_seats"
)

coalition_key <- function(S) {
  if (length(S) == 0) return("EMPTY")
  key <- paste(sort(S), collapse = "__")
  key <- gsub("[^A-Za-z0-9_\\-]", "_", key)
  key
}

coalition_to_params <- function(S) {
  pars <- character(0)
  
  if ("mode" %in% S)    pars <- c(pars, "asc_car", "asc_sm")
  if ("time" %in% S)    pars <- c(pars, "b_time")
  if ("cost" %in% S)    pars <- c(pars, "b_cost")
  if ("age" %in% S)     pars <- c(pars, "b_age")
  if ("ga" %in% S)      pars <- c(pars, "b_ga")
  if ("luggage" %in% S) pars <- c(pars, "b_luggage")
  if ("seats" %in% S)   pars <- c(pars, "b_seats")
  
  pars
}

fill_beta <- function(active_par, par_values) {
  beta <- setNames(rep(0, length(all_param_names)), all_param_names)
  if (length(active_par) > 0) beta[active_par] <- par_values
  beta
}

compute_utilities <- function(df, beta) {
  V_train <-
    beta["b_time"] * df$train_tt +
    beta["b_cost"] * df$train_co +
    beta["b_age"]  * df$age +
    beta["b_ga"]   * df$ga
  
  V_sm <-
    beta["asc_sm"] +
    beta["b_time"]  * df$sm_tt +
    beta["b_cost"]  * df$sm_co +
    beta["b_ga"]    * df$ga +
    beta["b_seats"] * df$sm_seats
  
  V_car <-
    beta["asc_car"] +
    beta["b_time"]    * df$car_tt +
    beta["b_cost"]    * df$car_co +
    beta["b_luggage"] * df$luggage
  
  cbind(train = V_train, sm = V_sm, car = V_car)
}

softmax_rows <- function(V) {
  maxV <- apply(V, 1, max)
  Vc <- V - maxV
  expV <- exp(Vc)
  denom <- rowSums(expV)
  expV / denom
}

loglik_mnl <- function(choice, probs) {
  idx <- cbind(seq_along(choice), choice)
  p <- probs[idx]
  p <- pmax(p, 1e-300)
  sum(log(p))
}

generate_all_coalitions <- function(players) {
  out <- list(character(0))
  if (length(players) == 0) return(out)
  
  for (k in seq_along(players)) {
    cmb <- combn(players, k, simplify = FALSE)
    out <- c(out, cmb)
  }
  out
}

cleanup_apollo_files <- function(model_name, out_dir) {
  files <- list.files(
    out_dir,
    pattern = paste0("^", model_name),
    full.names = TRUE
  )
  if (length(files) > 0) unlink(files, force = TRUE)
}

# -----------------------------
# 4) Estimación por coalición
# -----------------------------
estimate_mnl_coalition_apollo <- function(S, df, ll_null, benchmark_start) {
  
  key <- coalition_key(S)
  cache_file <- file.path("output/cache_coaliciones", paste0("coal_", key, ".rds"))
  
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  
  active_par <- coalition_to_params(S)
  
  # Coalición vacía = modelo nulo
  if (length(active_par) == 0) {
    beta_full <- fill_beta(character(0), numeric(0))
    V <- compute_utilities(df, beta_full)
    probs <- softmax_rows(V)
    ll_manual <- loglik_mnl(df$choice, probs)
    
    out <- list(
      coalition = S,
      key = key,
      active_par = active_par,
      beta = beta_full,
      ll_final = ll_manual,
      rho2 = 1 - (ll_manual / ll_null),
      probs = probs,
      converged = TRUE,
      method = "manual_empty"
    )
    
    saveRDS(out, cache_file)
    return(out)
  }
  
  database <- df
  model_name <- paste0("smcoal_", key)
  out_dir <- "output/apollo_coalitions"
  
  apollo_control <- list(
    modelName       = model_name,
    modelDescr      = paste("Swissmetro coalition", key),
    indivID         = "id",
    panelData       = TRUE,
    outputDirectory = out_dir
  )
  
  # Warm start desde el benchmark validado
  apollo_beta <- benchmark_start[active_par]
  apollo_fixed <- c()
  
  apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
    
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    
    V <- list(train = 0, sm = 0, car = 0)
    
    if ("asc_sm" %in% names(apollo_beta))  V[["sm"]]  <- V[["sm"]]  + asc_sm
    if ("asc_car" %in% names(apollo_beta)) V[["car"]] <- V[["car"]] + asc_car
    
    if ("b_time" %in% names(apollo_beta)) {
      V[["train"]] <- V[["train"]] + b_time * train_tt
      V[["sm"]]    <- V[["sm"]]    + b_time * sm_tt
      V[["car"]]   <- V[["car"]]   + b_time * car_tt
    }
    
    if ("b_cost" %in% names(apollo_beta)) {
      V[["train"]] <- V[["train"]] + b_cost * train_co
      V[["sm"]]    <- V[["sm"]]    + b_cost * sm_co
      V[["car"]]   <- V[["car"]]   + b_cost * car_co
    }
    
    if ("b_age" %in% names(apollo_beta)) {
      V[["train"]] <- V[["train"]] + b_age * age
    }
    
    if ("b_ga" %in% names(apollo_beta)) {
      V[["train"]] <- V[["train"]] + b_ga * ga
      V[["sm"]]    <- V[["sm"]]    + b_ga * ga
    }
    
    if ("b_luggage" %in% names(apollo_beta)) {
      V[["car"]] <- V[["car"]] + b_luggage * luggage
    }
    
    if ("b_seats" %in% names(apollo_beta)) {
      V[["sm"]] <- V[["sm"]] + b_seats * sm_seats
    }
    
    P <- list()
    
    mnl_settings <- list(
      alternatives = c(train = 1, sm = 2, car = 3),
      avail        = list(train = train_av, sm = sm_av, car = car_av),
      choiceVar    = choice,
      V            = V
    )
    
    P[["model"]] <- apollo_mnl(mnl_settings, functionality)
    P <- apollo_panelProd(P, apollo_inputs, functionality)
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    
    return(P)
  }
  
  apollo_inputs <- NULL
  modelo <- NULL
  
  invisible(capture.output({
    apollo_inputs <- apollo_validateInputs()
  }))
  
  modelo <- tryCatch(
    {
      invisible(capture.output({
        apollo_estimate(
          apollo_beta,
          apollo_fixed,
          apollo_probabilities,
          apollo_inputs
        ) -> modelo_tmp
      }))
      modelo_tmp
    },
    error = function(e) NULL
  )
  
  if (is.null(modelo)) {
    stop(paste("Apollo falló en la coalición:", key))
  }
  
  beta_active <- modelo$estimate
  beta_full <- fill_beta(active_par, beta_active)
  
  # Valor oficial del ajuste usando las mismas probs que usaremos luego
  V <- compute_utilities(df, beta_full)
  probs <- softmax_rows(V)
  ll_manual <- loglik_mnl(df$choice, probs)
  
  out <- list(
    coalition = S,
    key = key,
    active_par = active_par,
    beta = beta_full,
    ll_final = ll_manual,
    rho2 = 1 - (ll_manual / ll_null),
    probs = probs,
    converged = TRUE,
    method = "Apollo_warmstart"
  )
  
  saveRDS(out, cache_file)
  cleanup_apollo_files(model_name, out_dir)
  
  out
}

# -----------------------------
# 5) Estimar todas las coaliciones
# -----------------------------
coalitions <- generate_all_coalitions(players)

cat("Número total de coaliciones:", length(coalitions), "\n")
cat("Debería ser 2^P =", 2^P, "\n\n")

results_list <- vector("list", length(coalitions))
names(results_list) <- vapply(coalitions, coalition_key, character(1))

for (i in seq_along(coalitions)) {
  S <- coalitions[[i]]
  key <- coalition_key(S)
  
  cat(sprintf("[%3d/%3d] Estimando coalición: %s\n",
              i, length(coalitions), key))
  
  results_list[[key]] <- estimate_mnl_coalition_apollo(
    S = S,
    df = sm,
    ll_null = ll_null,
    benchmark_start = benchmark_start
  )
}

# -----------------------------
# 6) Chequeo modelo completo
# -----------------------------
full_key <- coalition_key(players)
full_res <- results_list[[full_key]]

benchmark_probs <- softmax_rows(compute_utilities(sm, benchmark_start))
benchmark_ll <- loglik_mnl(sm$choice, benchmark_probs)
benchmark_rho2 <- 1 - benchmark_ll / ll_null

comparison_full <- tibble(
  parametro = names(benchmark_start),
  benchmark = as.numeric(benchmark_start),
  estimate  = as.numeric(full_res$beta[names(benchmark_start)]),
  diff      = estimate - benchmark
)

cat("\n=====================================================\n")
cat("Chequeo de coalición completa\n")
cat("=====================================================\n")
cat("LL benchmark validado =", round(benchmark_ll, 3), "\n")
cat("LL full actual        =", round(full_res$ll_final, 3), "\n")
cat("rho2 benchmark        =", round(benchmark_rho2, 6), "\n")
cat("rho2 full actual      =", round(full_res$rho2, 6), "\n\n")
print(comparison_full)

write_csv(comparison_full, "output/tablas/swissmetro_full_model_vs_benchmark.csv")

# -----------------------------
# 7) Resumen de coaliciones
# -----------------------------
coalition_summary <- bind_rows(lapply(results_list, function(x) {
  tibble(
    key = x$key,
    coalition_size = length(x$coalition),
    coalition = if (length(x$coalition) == 0) "EMPTY" else paste(x$coalition, collapse = ", "),
    ll_final = x$ll_final,
    rho2 = x$rho2,
    converged = x$converged,
    method = x$method
  )
}))

write_csv(coalition_summary, "output/tablas/swissmetro_coalitions_summary.csv")

# -----------------------------
# 8) Shapley global exacto
# -----------------------------
shapley_global <- setNames(rep(0, P), players)

for (p in players) {
  absent_sets <- coalitions[!vapply(coalitions, function(S) p %in% S, logical(1))]
  
  for (S in absent_sets) {
    S_plus <- sort(c(S, p))
    
    key_S <- coalition_key(S)
    key_Sp <- coalition_key(S_plus)
    
    v_S <- results_list[[key_S]]$rho2
    v_Sp <- results_list[[key_Sp]]$rho2
    
    s <- length(S)
    weight <- factorial(s) * factorial(P - s - 1) / factorial(P)
    
    shapley_global[p] <- shapley_global[p] + weight * (v_Sp - v_S)
  }
}

shapley_global_df <- tibble(
  player = players,
  shapley_value = as.numeric(shapley_global),
  share_percent = 100 * shapley_value / sum(shapley_global)
)

cat("\n=====================================================\n")
cat("Shapley global exacto\n")
cat("=====================================================\n")
print(shapley_global_df)
cat("\n")
cat("Suma de Shapley global =", round(sum(shapley_global), 6), "\n")
cat("rho2 del modelo completo =", round(full_res$rho2, 6), "\n\n")

write_csv(shapley_global_df, "output/tablas/swissmetro_shapley_global_exact.csv")

# -----------------------------
# 9) SHAP local exacto
# -----------------------------
if (CALC_LOCAL) {
  
  shapley_local <- array(
    0,
    dim = c(N, 3, P),
    dimnames = list(NULL, c("train", "sm", "car"), players)
  )
  
  for (ip in seq_along(players)) {
    p <- players[ip]
    cat("Calculando SHAP local para jugador:", p, "\n")
    
    absent_sets <- coalitions[!vapply(coalitions, function(S) p %in% S, logical(1))]
    
    for (S in absent_sets) {
      S_plus <- sort(c(S, p))
      
      key_S <- coalition_key(S)
      key_Sp <- coalition_key(S_plus)
      
      probs_S <- results_list[[key_S]]$probs
      probs_Sp <- results_list[[key_Sp]]$probs
      
      s <- length(S)
      weight <- factorial(s) * factorial(P - s - 1) / factorial(P)
      
      shapley_local[, , ip] <- shapley_local[, , ip] + weight * (probs_Sp - probs_S)
    }
  }
  
  saveRDS(shapley_local, "output/swissmetro_shapley_local_exact.rds")
  
  probs_empty <- results_list[["EMPTY"]]$probs
  probs_full  <- results_list[[full_key]]$probs
  
  local_sum   <- apply(shapley_local, c(1, 2), sum)
  diff_check  <- local_sum - (probs_full - probs_empty)
  
  cat("\nMáximo error absoluto en eficiencia local:\n")
  print(max(abs(diff_check)))
  cat("\n")
}

# -----------------------------
# 10) Datos para gráficos
# -----------------------------
make_plot_df <- function(df, shapley_local) {
  bind_rows(
    tibble(alternative = "car",   feature = "CAR_CO",    feature_value = df$car_co,   shap = shapley_local[, "car",   "cost"]),
    tibble(alternative = "car",   feature = "CAR_TT",    feature_value = df$car_tt,   shap = shapley_local[, "car",   "time"]),
    tibble(alternative = "car",   feature = "LUGGAGE",   feature_value = df$luggage,  shap = shapley_local[, "car",   "luggage"]),
    
    tibble(alternative = "train", feature = "GA",        feature_value = df$ga,       shap = shapley_local[, "train", "ga"]),
    tibble(alternative = "train", feature = "TRAIN_TT",  feature_value = df$train_tt, shap = shapley_local[, "train", "time"]),
    tibble(alternative = "train", feature = "TRAIN_CO",  feature_value = df$train_co, shap = shapley_local[, "train", "cost"]),
    tibble(alternative = "train", feature = "AGE",       feature_value = df$age,      shap = shapley_local[, "train", "age"]),
    
    tibble(alternative = "sm",    feature = "SM_TT",     feature_value = df$sm_tt,    shap = shapley_local[, "sm",    "time"]),
    tibble(alternative = "sm",    feature = "GA",        feature_value = df$ga,       shap = shapley_local[, "sm",    "ga"]),
    tibble(alternative = "sm",    feature = "SM_CO",     feature_value = df$sm_co,    shap = shapley_local[, "sm",    "cost"]),
    tibble(alternative = "sm",    feature = "SM_SEATS",  feature_value = df$sm_seats, shap = shapley_local[, "sm",    "seats"])
  ) |>
    group_by(alternative, feature) |>
    mutate(mean_abs = mean(abs(shap), na.rm = TRUE)) |>
    ungroup()
}

plot_shap_summary_alt <- function(plot_df, alt_name, title_text = NULL) {
  
  dat <- plot_df |>
    filter(alternative == alt_name) |>
    group_by(feature) |>
    mutate(ord = mean(abs(shap), na.rm = TRUE)) |>
    ungroup()
  
  rng <- range(dat$feature_value, na.rm = TRUE)
  if (diff(rng) == 0) {
    dat <- dat |> mutate(color_value = 0.5)
  } else {
    dat <- dat |> mutate(color_value = (feature_value - rng[1]) / diff(rng))
  }
  
  order_levels <- dat |>
    distinct(feature, ord) |>
    arrange(ord) |>
    pull(feature)
  
  dat$feature <- factor(dat$feature, levels = order_levels)
  
  ggplot(dat, aes(x = shap, y = feature, color = color_value)) +
    ggbeeswarm::geom_quasirandom(
      width = 0.28,
      alpha = 0.75,
      size = 0.7,
      groupOnX = FALSE
    ) +
    scale_color_gradientn(
      colours = c("#2b8cbe", "#7bccc4", "#f03b87"),
      values = c(0, 0.5, 1),
      limits = c(0, 1),
      labels = c("Low", "High"),
      breaks = c(0, 1),
      name = "Feature value"
    ) +
    labs(
      title = ifelse(is.null(title_text), alt_name, title_text),
      x = "SHAP value (impact on model output)",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank()
    )
}

make_force_data <- function(row_id, alt_name, df, shapley_local, probs_empty) {
  
  if (alt_name == "car") {
    out <- tibble(
      feature = c("CAR_CO", "CAR_TT", "LUGGAGE"),
      value   = c(df$car_co[row_id], df$car_tt[row_id], df$luggage[row_id]),
      shap    = c(
        shapley_local[row_id, "car", "cost"],
        shapley_local[row_id, "car", "time"],
        shapley_local[row_id, "car", "luggage"]
      )
    )
    base0 <- probs_empty[row_id, "car"]
  }
  
  if (alt_name == "train") {
    out <- tibble(
      feature = c("GA", "TRAIN_TT", "TRAIN_CO", "AGE"),
      value   = c(df$ga[row_id], df$train_tt[row_id], df$train_co[row_id], df$age[row_id]),
      shap    = c(
        shapley_local[row_id, "train", "ga"],
        shapley_local[row_id, "train", "time"],
        shapley_local[row_id, "train", "cost"],
        shapley_local[row_id, "train", "age"]
      )
    )
    base0 <- probs_empty[row_id, "train"]
  }
  
  if (alt_name == "sm") {
    out <- tibble(
      feature = c("SM_TT", "GA", "SM_CO", "SM_SEATS"),
      value   = c(df$sm_tt[row_id], df$ga[row_id], df$sm_co[row_id], df$sm_seats[row_id]),
      shap    = c(
        shapley_local[row_id, "sm", "time"],
        shapley_local[row_id, "sm", "ga"],
        shapley_local[row_id, "sm", "cost"],
        shapley_local[row_id, "sm", "seats"]
      )
    )
    base0 <- probs_empty[row_id, "sm"]
  }
  
  out <- out |>
    mutate(
      abs_shap = abs(shap),
      sign = ifelse(shap >= 0, "positive", "negative")
    ) |>
    arrange(desc(abs_shap))
  
  end_vec <- base0 + cumsum(out$shap)
  start_vec <- c(base0, head(end_vec, -1))
  
  out |>
    mutate(
      baseline = base0,
      start = start_vec,
      end = end_vec,
      mid = (start + end) / 2
    )
}

plot_force_like_alt <- function(row_id, alt_name, df, shapley_local, probs_empty) {
  
  dat <- make_force_data(row_id, alt_name, df, shapley_local, probs_empty)
  baseline <- dat$baseline[1]
  pred <- baseline + sum(dat$shap)
  
  ggplot(dat) +
    geom_segment(
      aes(x = start, xend = end, y = feature, yend = feature, color = sign),
      linewidth = 8,
      lineend = "butt"
    ) +
    geom_vline(xintercept = baseline, linetype = "dashed", color = "grey40") +
    geom_text(
      aes(x = mid, y = feature, label = paste0(feature, " = ", round(value, 3))),
      color = "white",
      size = 3
    ) +
    scale_color_manual(values = c(negative = "#2b8cbe", positive = "#f03b87")) +
    labs(
      title = paste0(toupper(alt_name), " | base = ", round(baseline, 3), " | pred = ", round(pred, 3)),
      x = "Contribución acumulada",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )
}

# -----------------------------
# 11) Crear gráficos SHAP
# -----------------------------
if (CALC_LOCAL && MAKE_PLOTS) {
  
  plot_df <- make_plot_df(sm, shapley_local)
  
  p_car <- plot_shap_summary_alt(plot_df, "car",   "Car")
  p_sm  <- plot_shap_summary_alt(plot_df, "sm",    "Swissmetro")
  p_tr  <- plot_shap_summary_alt(plot_df, "train", "Train")
  
  p_summary_all <- (p_car | p_sm) / p_tr +
    plot_annotation(title = "SHAP summary plots by alternative")
  
  ggsave(
    filename = "output/figuras/swissmetro_shap_summary_by_alt.png",
    plot = p_summary_all,
    width = 12,
    height = 9,
    dpi = 300
  )
  
  row_id_force <- 1
  
  p_force_car <- plot_force_like_alt(row_id_force, "car", sm, shapley_local, probs_empty)
  p_force_sm  <- plot_force_like_alt(row_id_force, "sm", sm, shapley_local, probs_empty)
  p_force_tr  <- plot_force_like_alt(row_id_force, "train", sm, shapley_local, probs_empty)
  
  p_force_all <- p_force_car / p_force_sm / p_force_tr +
    plot_annotation(title = paste0("Local explanation for row ", row_id_force))
  
  ggsave(
    filename = "output/figuras/swissmetro_force_like_row1.png",
    plot = p_force_all,
    width = 12,
    height = 10,
    dpi = 300
  )
}

# -----------------------------
# 12) Curva de predicción con Apollo
# -----------------------------
# Se usa apollo_prediction() y, si su salida no es utilizable
# en esta sesión, se recurre al cálculo manual con softmax.

database <- sm

apollo_control <- list(
  modelName       = "Swissmetro_prediction_plot",
  modelDescr      = "Prediction plot for Swissmetro",
  indivID         = "id",
  panelData       = TRUE,
  outputDirectory = "output/apollo_coalitions"
)

apollo_beta  <- benchmark_start
apollo_fixed <- c()
apollo_inputs <- apollo_validateInputs()

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
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
  
  P <- list()
  
  mnl_settings <- list(
    alternatives = c(train = 1, sm = 2, car = 3),
    avail        = list(train = train_av, sm = sm_av, car = car_av),
    choiceVar    = choice,
    V            = V
  )
  
  P[["model"]] <- apollo_mnl(mnl_settings, functionality)
  P <- apollo_panelProd(P, apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  
  return(P)
}

extract_prediction_matrix <- function(pred_obj) {
  
  if (is.matrix(pred_obj) || is.data.frame(pred_obj)) {
    pred_df <- as.data.frame(pred_obj)
    nm <- tolower(names(pred_df))
    if (all(c("train", "sm", "car") %in% nm)) {
      names(pred_df) <- nm
      return(pred_df[, c("train", "sm", "car")])
    }
  }
  
  if (is.list(pred_obj)) {
    # caso: lista con componente "model"
    if ("model" %in% names(pred_obj)) {
      return(extract_prediction_matrix(pred_obj$model))
    }
    
    # caso: lista nombrada con train/sm/car
    nm <- tolower(names(pred_obj))
    if (all(c("train", "sm", "car") %in% nm)) {
      pred_df <- as.data.frame(pred_obj)
      names(pred_df) <- nm
      return(pred_df[, c("train", "sm", "car")])
    }
  }
  
  NULL
}

sm_original <- sm

tt_grid <- seq(
  from = quantile(sm$sm_tt, 0.05, na.rm = TRUE),
  to   = quantile(sm$sm_tt, 0.95, na.rm = TRUE),
  length.out = 30
)

pred_curve <- lapply(tt_grid, function(tt_val) {
  
  database <<- sm_original
  database$sm_tt <- tt_val
  
  apollo_inputs_tmp <- apollo_validateInputs()
  
  preds_obj <- tryCatch(
    apollo_prediction(
      benchmark_start,
      apollo_probabilities,
      apollo_inputs_tmp,
      prediction_settings = list(summary = FALSE)
    ),
    error = function(e) NULL
  )
  
  preds_df <- extract_prediction_matrix(preds_obj)
  
  # Fallback manual si Apollo no devuelve algo usable
  if (is.null(preds_df)) {
    V_tmp <- compute_utilities(database, benchmark_start)
    preds_df <- as.data.frame(softmax_rows(V_tmp))
    names(preds_df) <- c("train", "sm", "car")
  }
  
  tibble(
    sm_tt = tt_val,
    mean_prob_train = mean(preds_df$train),
    mean_prob_sm    = mean(preds_df$sm),
    mean_prob_car   = mean(preds_df$car)
  )
}) |>
  bind_rows()

pred_curve_long <- pred_curve |>
  pivot_longer(
    cols = starts_with("mean_prob_"),
    names_to = "alternative",
    values_to = "probability"
  ) |>
  mutate(
    alternative = dplyr::recode(
      alternative,
      mean_prob_train = "Train",
      mean_prob_sm    = "Swissmetro",
      mean_prob_car   = "Car"
    )
  )

p_pred <- ggplot(pred_curve_long, aes(x = sm_tt, y = probability, color = alternative)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Predicted probabilities as SM_TT varies",
    x = "Swissmetro travel time",
    y = "Average predicted probability",
    color = "Alternative"
  ) +
  theme_minimal(base_size = 11)

ggsave(
  filename = "output/figuras/swissmetro_prediction_curve_sm_tt.png",
  plot = p_pred,
  width = 9,
  height = 5,
  dpi = 300
)

database <- sm_original

# -----------------------------
# 13) Guardar objetos finales
# -----------------------------
final_object <- list(
  players = players,
  ll_null = ll_null,
  benchmark_start = benchmark_start,
  coalition_summary = coalition_summary,
  full_model = full_res,
  shapley_global = shapley_global_df,
  shapley_local = if (CALC_LOCAL) shapley_local else NULL
)

saveRDS(final_object, "output/swissmetro_shapley_exact_results.rds")

cat("=====================================================\n")
cat("Proceso terminado\n")
cat("=====================================================\n")
cat("Archivos clave:\n")
cat("- output/tablas/swissmetro_coalitions_summary.csv\n")
cat("- output/tablas/swissmetro_full_model_vs_benchmark.csv\n")
cat("- output/tablas/swissmetro_shapley_global_exact.csv\n")
if (CALC_LOCAL) {
  cat("- output/swissmetro_shapley_local_exact.rds\n")
  cat("- output/figuras/swissmetro_shap_summary_by_alt.png\n")
  cat("- output/figuras/swissmetro_force_like_row1.png\n")
}
cat("- output/figuras/swissmetro_prediction_curve_sm_tt.png\n")
cat("- output/swissmetro_shapley_exact_results.rds\n")