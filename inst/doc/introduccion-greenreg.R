## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup--------------------------------------------------------------------
library(GREENREG)

## -----------------------------------------------------------------------------
data(datos_rendimiento_maiz)

modelo <- greenreg_rlm(rendimiento_maiz_ton_ha ~ ., data = datos_rendimiento_maiz)
modelo

## ----fig.alt = "Diagnosticos graficos del modelo de regresion lineal multiple"----
plot(modelo)

## -----------------------------------------------------------------------------
data(datos_roya_cafe)

modelo_log <- greenreg_logistico(
  presencia_roya ~ humedad_relativa_pct + altitud_msnm + manejo_sombra,
  data = datos_roya_cafe
)
modelo_log

## -----------------------------------------------------------------------------
data(datos_anomalia_temperatura)

modelo_ts <- greenreg_modelo_ar(datos_anomalia_temperatura$anomalia_c, p = 2)
modelo_ts

## -----------------------------------------------------------------------------
greenreg_analisis_datos(datos_rendimiento_maiz, "rendimiento_maiz_ton_ha")

