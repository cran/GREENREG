# Pruebas de estres y casos limite adaptadas de 02_estres.R. Un paquete
# publicado debe fallar con un mensaje claro o manejar el caso; nunca abortar
# con un error interno de R. Estas pruebas son la regresion de C-02, C-03,
# C-04 y C-05 (todas citadas como ERROR en el informe antes del fix).

silencio <- function(expr) invisible(capture.output(suppressMessages(suppressWarnings(expr))))

test_that("C-02: NA en los predictores no rompe data_plot (greenreg_rls, greenreg_rlm, plot)", {
  d_na <- datos_rendimiento_maiz
  d_na$precipitacion_mm[c(3, 10)] <- NA
  d_na$temperatura_c[7] <- NA

  expect_no_error(silencio(print(greenreg_rls(rendimiento_maiz_ton_ha ~ precipitacion_mm, d_na))))
  expect_no_error(silencio(print(greenreg_rlm(rendimiento_maiz_ton_ha ~ ., d_na))))
  expect_no_error({
    pdf(NULL); on.exit(dev.off(), add = TRUE)
    silencio(plot(greenreg_rlm(rendimiento_maiz_ton_ha ~ ., d_na)))
  })
})

test_that("C-02: el aviso de filas omitidas por NA aparece cuando corresponde", {
  d_na <- datos_rendimiento_maiz
  d_na$precipitacion_mm[c(3, 10)] <- NA
  s <- paste(capture.output(print(greenreg_rlm(rendimiento_maiz_ton_ha ~ ., d_na))), collapse = "\n")
  expect_true(grepl("omitieron.*NA", s, ignore.case = TRUE))
})

test_that("C-03: greenreg_logistico() con un solo predictor no aborta (VIF no aplicable)", {
  expect_no_error(silencio(print(
    greenreg_logistico(presencia_roya ~ humedad_relativa_pct, data = datos_roya_cafe)
  )))
})

test_that("C-03: separacion completa no aborta greenreg_logistico()", {
  d_sep <- data.frame(y = rep(0:1, each = 30), x = c(rnorm(30, 0), rnorm(30, 10)))
  expect_no_error(silencio(print(greenreg_logistico(y ~ x, d_sep))))
})

test_that("C-04: colinealidad perfecta no rompe greenreg_rlm() ni su plot()", {
  d_col <- datos_rendimiento_maiz
  d_col$precip2 <- d_col$precipitacion_mm * 2  # colinealidad perfecta

  expect_no_error(silencio(print(greenreg_rlm(rendimiento_maiz_ton_ha ~ ., d_col))))
  expect_no_error({
    pdf(NULL); on.exit(dev.off(), add = TRUE)
    silencio(plot(greenreg_rlm(rendimiento_maiz_ton_ha ~ ., d_col)))
  })
})

test_that("C-04: la colinealidad perfecta se reporta como alerta explicita", {
  d_col <- datos_rendimiento_maiz
  d_col$precip2 <- d_col$precipitacion_mm * 2
  s <- paste(capture.output(print(greenreg_rlm(rendimiento_maiz_ton_ha ~ ., d_col))), collapse = "\n")
  expect_true(grepl("PELIGRO.*[Cc]olinealidad", s))
})

test_that("VIF alto (pero no perfecto) sigue generando la alerta esperada", {
  d_col2 <- datos_rendimiento_maiz
  d_col2$precip_ruido <- d_col2$precipitacion_mm + rnorm(nrow(d_col2), 0, 1)
  s <- paste(capture.output(print(greenreg_rlm(rendimiento_maiz_ton_ha ~ ., d_col2))), collapse = "\n")
  expect_true(grepl("PELIGRO|VIF > 10", s))
})

test_that("C-05: n pequenio no rompe greenreg_rlm() ni greenreg_modelo_ar()", {
  expect_no_error(silencio(print(greenreg_rlm(rendimiento_maiz_ton_ha ~ ., datos_rendimiento_maiz[1:5, ]))))
  expect_no_error(silencio(print(greenreg_rlm(rendimiento_maiz_ton_ha ~ ., datos_rendimiento_maiz[1:4, ]))))
  set.seed(2026)
  expect_no_error(silencio(print(greenreg_modelo_ar(rnorm(5), p = 2))))
})

test_that("greenreg_modelo_ma con q = 0 no aborta", {
  set.seed(2026)
  expect_no_error(silencio(print(greenreg_modelo_ma(rnorm(50), q = 0))))
})

test_that("greenreg_rls con mas de un predictor da un error claro que remite a greenreg_rlm (C-14)", {
  expect_error(
    silencio(greenreg_rls(rendimiento_maiz_ton_ha ~ precipitacion_mm + temperatura_c, datos_rendimiento_maiz)),
    "un solo predictor"
  )
})
