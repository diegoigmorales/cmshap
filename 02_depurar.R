# =========================================================
# 02_depurar.R
# Limpieza de Swissmetro usando solo R Base
# =========================================================

# 1. Cargar swissmetro.dat
archivo_entrada <- "swissmetro.dat"
df <- read.table(
  archivo_entrada,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# 2. Eliminar observaciones con CHOICE faltante o invalido
df <- df[!is.na(df$CHOICE) & df$CHOICE %in% c(1, 2, 3), ]

# 3. Eliminar observaciones sin disponibilidad simultanea
df <- df[df$CAR_AV == 1 & df$TRAIN_AV == 1 & df$SM_AV == 1, ]

# 4. Conservar solo variables del estudio:
# CHOICE, TT, COST, AGE, GA, SEATS y LUGGAGE
df <- data.frame(
  CHOICE = df$CHOICE,
  CAR_TT = df$CAR_TT,
  TRAIN_TT = df$TRAIN_TT,
  SM_TT = df$SM_TT,
  CAR_CO = df$CAR_CO,
  TRAIN_CO = df$TRAIN_CO,
  SM_CO = df$SM_CO,
  AGE = df$AGE,
  GA = df$GA,
  SM_SEATS = df$SM_SEATS,
  LUGGAGE = df$LUGGAGE
)

# 5. Exportar resultado final
write.csv(df, "swissmetro_limpio.csv", row.names = FALSE)
