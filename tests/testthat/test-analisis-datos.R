# Pruebas de greenreg_analisis_datos() (clasificacion de la variable respuesta),
# adaptadas de los checks ad-01 a ad-08. ad-04/ad-05/ad-08 (hallazgos C-13/C-12/C-11)
# estan activas: verifican las correcciones.

test_that("greenreg_analisis_datos: clasifica correctamente continua/binaria/conteo", {
  r_cont <- capture.output(suppressMessages(
    greenreg_analisis_datos(datos_rendimiento_maiz, "rendimiento_maiz_ton_ha")))
  r_bin <- capture.output(suppressMessages(
    greenreg_analisis_datos(datos_roya_cafe, "presencia_roya")))
  r_cnt <- capture.output(suppressMessages(
    greenreg_analisis_datos(datos_tala_ilegal_poisson, "numero_talas")))

  expect_true(any(grepl("Sugerencia: Opción 1", r_cont)))
  expect_true(any(grepl("Sugerencia: Opción 2", r_bin)))
  expect_true(any(grepl("Sugerencia: Opción 2", r_cnt)))
})

test_that("greenreg_analisis_datos: detecta serie de tiempo con columna de fecha en texto", {
  r <- capture.output(suppressMessages(
    greenreg_analisis_datos(datos_anomalia_temperatura, "anomalia_c")))
  expect_true(any(grepl("Sugerencia: Opción 3", r)))
})

test_that("greenreg_analisis_datos: no clasifica como serie una variable ruido blanco con fecha", {
  set.seed(1)
  d <- data.frame(fecha = as.character(seq(as.Date("2000-01-01"), by = "month", length.out = 120)),
                   y = rnorm(120))
  r <- capture.output(suppressMessages(greenreg_analisis_datos(d, "y")))
  expect_false(any(grepl("Sugerencia: Opción 3", r)))
})

test_that("greenreg_analisis_datos: Opcion B muestra la formula seleccionada, no siempre '~ .' (C-13)", {
  r <- capture.output(suppressMessages(
    greenreg_analisis_datos(datos_rendimiento_maiz, "rendimiento_maiz_ton_ha")))
  expect_false(any(grepl("modelo_rlm_opt <- greenreg_rlm\\(.* ~ \\.,", r)))
})

test_that("greenreg_analisis_datos: detecta serie de tiempo cuando la columna es de clase ts (C-12)", {
  set.seed(12)
  d <- data.frame(t = 1:60)
  d$y <- ts(cumsum(rnorm(60)))   # caminata aleatoria: autocorrelacion fuerte, columna de clase ts
  r <- capture.output(suppressMessages(greenreg_analisis_datos(d, "y")))
  expect_true(any(grepl("Sugerencia: Opción 3", r)))
})

test_that("greenreg_analisis_datos: no clasifica como conteo (Poisson) una variable continua con valores enteros (C-11)", {
  # El menu siempre menciona "Poisson": se revisa la SUGERENCIA, no todo el texto.
  set.seed(11)
  d <- data.frame(y = round(rnorm(100, 170, 15)), x = rnorm(100))
  r <- capture.output(suppressMessages(greenreg_analisis_datos(d, "y")))
  expect_true(any(grepl("Sugerencia: Opción 1", r)))
  expect_false(any(grepl("Sugerencia: Opción 2", r)))
})
