# Pruebas de los alias obsoletos (E-03): rlm(), rls(), logistico(),
# reg_poisson(), modelo_ar/ma/arma/arima() y analisis_datos() deben seguir
# funcionando exactamente igual que sus versiones greenreg_*, pero avisando
# con .Deprecated() (rlm() en particular enmascaraba MASS::rlm() en
# silencio; ahora al menos avisa).

test_that("rlm() (obsoleta) da el mismo resultado que greenreg_rlm() y avisa", {
  expect_warning(
    m_viejo <- rlm(rendimiento_maiz_ton_ha ~ ., data = datos_rendimiento_maiz),
    class = "deprecatedWarning"
  )
  m_nuevo <- greenreg_rlm(rendimiento_maiz_ton_ha ~ ., data = datos_rendimiento_maiz)
  expect_equal(coef(m_viejo$modelo), coef(m_nuevo$modelo))
})

test_that("rls() (obsoleta) da el mismo resultado que greenreg_rls() y avisa", {
  expect_warning(
    m_viejo <- rls(rendimiento_maiz_ton_ha ~ precipitacion_mm, data = datos_rendimiento_maiz),
    class = "deprecatedWarning"
  )
  m_nuevo <- greenreg_rls(rendimiento_maiz_ton_ha ~ precipitacion_mm, data = datos_rendimiento_maiz)
  expect_equal(coef(m_viejo$modelo), coef(m_nuevo$modelo))
})

test_that("modelo_ar() (obsoleta) da el mismo resultado que greenreg_modelo_ar() y avisa", {
  serie <- datos_anomalia_temperatura$anomalia_c
  expect_warning(m_viejo <- modelo_ar(serie, p = 2), class = "deprecatedWarning")
  m_nuevo <- greenreg_modelo_ar(serie, p = 2)
  expect_equal(coef(m_viejo$modelo), coef(m_nuevo$modelo))
})

test_that("analisis_datos() (obsoleta) llama a greenreg_analisis_datos()", {
  expect_warning(
    analisis_datos(datos_rendimiento_maiz, "rendimiento_maiz_ton_ha", verbose = FALSE),
    class = "deprecatedWarning"
  )
})
