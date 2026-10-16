# Pruebas de greenreg_rls() (regresion lineal simple) contra stats::lm() y lmtest,
# adaptadas de los checks greenreg_rls-01 a greenreg_rls-10 del informe de revision de codigo
# (Israel Lerma Serna, DICIFO-UACh, septiembre 2026).

tol <- 1e-8

salida <- function(obj) paste(utils::capture.output(print(obj)), collapse = "\n")

test_that("greenreg_rls: coeficientes iguales a lm()", {
  m  <- greenreg_rls(rendimiento_maiz_ton_ha ~ precipitacion_mm, data = datos_rendimiento_maiz)
  mr <- stats::lm(rendimiento_maiz_ton_ha ~ precipitacion_mm, data = datos_rendimiento_maiz)
  expect_true(all(abs(coef(m$modelo) - coef(mr)) < tol))
})

test_that("greenreg_rls: R2, Shapiro, Breusch-Pagan y Durbin-Watson impresos coinciden con las referencias", {
  m  <- greenreg_rls(rendimiento_maiz_ton_ha ~ precipitacion_mm, data = datos_rendimiento_maiz)
  mr <- stats::lm(rendimiento_maiz_ton_ha ~ precipitacion_mm, data = datos_rendimiento_maiz)
  s  <- salida(m)

  expect_true(grepl(sprintf("%.2f", summary(mr)$r.squared * 100), s, fixed = TRUE))
  expect_true(grepl(sprintf("%.4f", stats::shapiro.test(residuals(mr))$p.value), s, fixed = TRUE))
  expect_true(grepl(sprintf("%.4f", lmtest::bptest(mr)$p.value), s, fixed = TRUE))
  expect_true(grepl(sprintf("%.2f", lmtest::dwtest(mr)$statistic), s, fixed = TRUE))
})

test_that("greenreg_rls: expone el objeto lm subyacente", {
  m <- greenreg_rls(rendimiento_maiz_ton_ha ~ precipitacion_mm, data = datos_rendimiento_maiz)
  expect_s3_class(m$modelo, "lm")
})

# greenreg_rls-07/greenreg_rls-09 (informe) / E-12: decision de diseno. El paquete usa
# "OKEY" / "ALERTA" como veredictos y redacta "no se rechaza ..." en vez de afirmar
# la hipotesis nula ("parecen normales").
test_that("greenreg_rls: veredictos OKEY/ALERTA y redaccion 'no se rechaza' sin afirmar la hipotesis nula (E-12)", {
  s <- salida(greenreg_rls(rendimiento_maiz_ton_ha ~ precipitacion_mm, data = datos_rendimiento_maiz))
  expect_true(grepl("OKEY", s))
  expect_false(grepl("(^|[^A-Z])OK([^A-Z]|$)", s))   # sin "OK" suelto
  expect_false(grepl("parecen normales", s))
  expect_true(grepl("No se rechaza", s))
})

# greenreg_rls-08/E-07 (informe): 0% de las cadenas del paquete llevan acentos. Hallazgo
# E-07 (Alta): las cadenas de salida llevan acentos.
test_that("greenreg_rls: la salida incluye acentos correctos (E-07)", {
  s <- salida(greenreg_rls(rendimiento_maiz_ton_ha ~ precipitacion_mm, data = datos_rendimiento_maiz))
  expect_true(grepl("REGRESI\u00d3N|Regresi\u00f3n", s))
})
