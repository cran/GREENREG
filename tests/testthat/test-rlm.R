# Pruebas de greenreg_rlm() (regresion lineal multiple) contra stats::lm(), car::vif() y
# lmtest, adaptadas de los checks greenreg_rlm-01 a greenreg_rlm-06 del informe de revision de
# codigo. greenreg_rlm-05/greenreg_rlm-06 son la prueba de regresion para el hallazgo C-01.

tol <- 1e-6

salida <- function(obj) paste(utils::capture.output(print(obj)), collapse = "\n")

test_that("greenreg_rlm: coeficientes iguales a lm()", {
  m  <- greenreg_rlm(rendimiento_maiz_ton_ha ~ ., data = datos_rendimiento_maiz)
  mr <- stats::lm(rendimiento_maiz_ton_ha ~ ., data = datos_rendimiento_maiz)
  expect_true(all(abs(coef(m$modelo) - coef(mr)) < tol))
})

test_that("greenreg_rlm: VIF manual coincide con car::vif()", {
  skip_if_not_installed("car")
  m  <- greenreg_rlm(rendimiento_maiz_ton_ha ~ ., data = datos_rendimiento_maiz)
  mr <- stats::lm(rendimiento_maiz_ton_ha ~ ., data = datos_rendimiento_maiz)
  vif_ref <- car::vif(mr)
  expect_true(all(abs(m$vif[names(vif_ref)] - vif_ref) < 1e-6))
})

test_that("greenreg_rlm: R2 ajustado y Rainbow test impresos coinciden con las referencias", {
  m  <- greenreg_rlm(rendimiento_maiz_ton_ha ~ ., data = datos_rendimiento_maiz)
  mr <- stats::lm(rendimiento_maiz_ton_ha ~ ., data = datos_rendimiento_maiz)
  s  <- salida(m)
  expect_true(grepl(sprintf("%.2f", summary(mr)$adj.r.squared * 100), s, fixed = TRUE))
  expect_true(grepl(sprintf("%.4f", lmtest::raintest(mr)$p.value), s, fixed = TRUE))
})

test_that("C-01: la linea de Linealidad (Rainbow Test) no imprime el mensaje del VIF", {
  m <- greenreg_rlm(rendimiento_maiz_ton_ha ~ ., data = datos_rendimiento_maiz)
  s <- salida(m)
  lin_line <- grep("Rainbow", strsplit(s, "\n")[[1]], value = TRUE)
  expect_false(grepl("VIF|Multicolinealidad", lin_line))
})

test_that("C-01: notas$supuestos tiene 5 entradas con nombre (norm, homo, indep, lineal, vif)", {
  m <- greenreg_rlm(rendimiento_maiz_ton_ha ~ ., data = datos_rendimiento_maiz)
  expect_length(m$notas$supuestos, 5)
  expect_setequal(names(m$notas$supuestos), c("norm", "homo", "indep", "lineal", "vif"))
})
