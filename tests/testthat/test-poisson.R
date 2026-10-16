# Pruebas de greenreg_reg_poisson() contra stats::glm(poisson) y AER, adaptadas de los
# checks poi-01 a poi-07. poi-04/05/06 (hallazgos E-11/C-07) estan activas:
# verifican las correcciones.

salida <- function(obj) paste(utils::capture.output(print(obj)), collapse = "\n")

modelo_poi <- function() {
  greenreg_reg_poisson(numero_talas ~ densidad_poblacion + superficie_forestal_ha + dias_vigilancia,
              data = datos_tala_ilegal_poisson)
}
glm_ref <- function() {
  glm(numero_talas ~ densidad_poblacion + superficie_forestal_ha + dias_vigilancia,
      data = datos_tala_ilegal_poisson, family = poisson)
}

test_that("greenreg_reg_poisson: coeficientes iguales a glm(poisson)", {
  m <- modelo_poi(); mr <- glm_ref()
  expect_true(all(abs(coef(m$modelo) - coef(mr)) < 1e-6))
})

test_that("greenreg_reg_poisson: IRR = exp(coef)", {
  m <- modelo_poi(); mr <- glm_ref()
  expect_true(all(abs(m$irr_table[, "IRR"] - exp(coef(mr))) < 1e-6))
})

test_that("greenreg_reg_poisson: dispersion = Pearson chi2 / gl (lo que dice el codigo)", {
  m <- modelo_poi(); mr <- glm_ref()
  disp_pearson <- sum(residuals(mr, "pearson")^2) / df.residual(mr)
  expect_true(abs(m$dispersion - disp_pearson) < 1e-8)
})

test_that("greenreg_reg_poisson: la conclusion de dispersion aplica el umbral 1.5 sobre Pearson chi2/gl (E-11)", {
  m <- modelo_poi(); s <- salida(m)
  expect_equal(m$dispersion <= 1.5, grepl("OKEY:? Se cumple", s))
  expect_equal(m$dispersion > 1.5, grepl("SOBREDISPERSI", s))
})

test_that("greenreg_reg_poisson: la ecuacion impresa no trunca coeficientes a 0 y usa log(mu) (C-07)", {
  s <- salida(modelo_poi())
  expect_false(grepl("\\(0 \\*", s))
  expect_false(grepl("log\\(Y\\)", s))
})
