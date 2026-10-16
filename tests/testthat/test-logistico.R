# Pruebas de greenreg_logistico() (regresion logistica) contra stats::glm(binomial) y
# pROC, adaptadas de los checks log-01 a log-10 del informe de revision de
# codigo. log-05/06/08/10 (hallazgos C-08/C-09) estan activas: verifican las
# correcciones.

salida <- function(obj) paste(utils::capture.output(print(obj)), collapse = "\n")

modelo_log <- function() {
  greenreg_logistico(presencia_roya ~ humedad_relativa_pct + altitud_msnm + manejo_sombra,
            data = datos_roya_cafe)
}
glm_ref <- function() {
  glm(presencia_roya ~ humedad_relativa_pct + altitud_msnm + manejo_sombra,
      data = datos_roya_cafe, family = binomial)
}

test_that("greenreg_logistico: coeficientes iguales a glm(binomial)", {
  m <- modelo_log(); mr <- glm_ref()
  expect_true(all(abs(coef(m$modelo) - coef(mr)) < 1e-6))
})

test_that("greenreg_logistico: OR = exp(coef)", {
  m <- modelo_log(); mr <- glm_ref()
  expect_true(all(abs(m$odds_ratios[, "OR"] - exp(coef(mr))) < 1e-6))
})

test_that("greenreg_logistico: McFadden impreso coincide con la formula", {
  m <- modelo_log(); mr <- glm_ref(); s <- salida(m)
  mcf <- 1 - as.numeric(logLik(mr)) / as.numeric(logLik(update(mr, . ~ 1)))
  expect_true(grepl(sprintf("%.3f", mcf), s, fixed = TRUE))
})

test_that("greenreg_logistico: matriz de confusion coincide", {
  m <- modelo_log(); mr <- glm_ref()
  tc <- table(Real = mr$y, Predicho = as.integer(fitted(mr) > 0.5))
  expect_true(all(m$confusion == tc))
})

test_that("greenreg_logistico: AUC impreso = pROC::auc", {
  skip_if_not_installed("pROC")
  m <- modelo_log(); mr <- glm_ref(); s <- salida(m)
  auc_ref <- as.numeric(pROC::auc(mr$y, fitted(mr), quiet = TRUE))
  expect_true(grepl(sprintf("%.2f", auc_ref), s, fixed = TRUE))
})

test_that("greenreg_logistico: LRT impreso = anova(nulo, modelo)", {
  m <- modelo_log(); mr <- glm_ref(); s <- salida(m)
  lrt <- anova(update(mr, . ~ 1), mr, test = "Chisq")
  expect_true(grepl(sprintf("%.3f", lrt$Deviance[2]), s, fixed = TRUE))
})

test_that("greenreg_logistico: reporta sensibilidad/especificidad y tasa de no-informacion (C-09)", {
  s <- salida(modelo_log())
  expect_true(grepl("Sensibilidad|Especificidad", s, ignore.case = TRUE))
  expect_true(grepl("no informaci|mayoritaria|base", s, ignore.case = TRUE))
})

test_that("greenreg_logistico: el OR del intercepto no se imprime como 0.0000 (C-08)", {
  s <- salida(modelo_log())
  expect_false(grepl("\\(Intercept\\)\\s+0\\.0000\\s+0\\.0000", s))
})

test_that("greenreg_logistico: la etiqueta 'excelente' del McFadden respeta el rango 0.2-0.4 (C-09)", {
  m <- modelo_log(); mr <- glm_ref(); s <- salida(m)
  mcf <- 1 - as.numeric(logLik(mr)) / as.numeric(logLik(update(mr, . ~ 1)))
  expect_false(mcf > 0.4 && grepl("excelente", s))
})

test_that("C-03: greenreg_logistico() no aborta con un solo predictor", {
  expect_no_error(
    print(greenreg_logistico(presencia_roya ~ humedad_relativa_pct, data = datos_roya_cafe))
  )
})
