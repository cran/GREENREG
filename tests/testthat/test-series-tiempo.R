# Pruebas de greenreg_modelo_ar()/greenreg_modelo_ma()/greenreg_modelo_arma()/greenreg_modelo_arima() contra
# stats::arima(), adaptadas de los checks ts-01 a ts-16. ts-13/ts-14 (hallazgo
# C-10) estan activas.

salida <- function(obj) paste(utils::capture.output(print(obj)), collapse = "\n")

x <- datos_anomalia_temperatura$anomalia_c
y <- datos_nivel_presa$nivel_m

test_that("greenreg_modelo_ar/ma/arma/arima: coeficientes iguales a stats::arima()", {
  a  <- greenreg_modelo_ar(x, p = 2);                    ar_ref  <- arima(x, order = c(2, 0, 0))
  b  <- greenreg_modelo_ma(x, q = 1, include_mean = TRUE); ma_ref  <- arima(x, order = c(0, 0, 1))
  cc <- greenreg_modelo_arma(x, p = 1, q = 1, include_mean = TRUE); arma_ref <- arima(x, order = c(1, 0, 1))
  d  <- greenreg_modelo_arima(y, p = 1, d = 1, q = 1);    arima_ref <- arima(y, order = c(1, 1, 1))

  expect_true(all(abs(coef(a$modelo)  - coef(ar_ref))   < 1e-6))
  expect_true(all(abs(coef(b$modelo)  - coef(ma_ref))   < 1e-6))
  expect_true(all(abs(coef(cc$modelo) - coef(arma_ref)) < 1e-6))
  expect_true(all(abs(coef(d$modelo)  - coef(arima_ref)) < 1e-6))
})

test_that("greenreg_modelo_ar: modulos de raices inversas impresos y coinciden con forecast::arroots", {
  skip_if_not_installed("forecast")
  a <- greenreg_modelo_ar(x, p = 2); ar_ref <- arima(x, order = c(2, 0, 0)); sa <- salida(a)
  inv <- Mod(1 / polyroot(c(1, -coef(ar_ref)[1:2])))
  expect_true(all(sapply(sprintf("%.3f", inv), grepl, sa, fixed = TRUE)))

  fr <- forecast::Arima(x, order = c(2, 0, 0))
  inv_fr <- Mod(1 / polyroot(c(1, -coef(fr)[1:2])))
  expect_true(all(abs(sort(inv_fr) - sort(inv)) < 1e-6))
})

test_that("C-06: Ljung-Box usa fitdf = numero de parametros ARMA estimados", {
  a <- greenreg_modelo_ar(x, p = 2); ar_ref <- arima(x, order = c(2, 0, 0)); sa <- salida(a)

  # Reproduce exactamente la formula de lag/fitdf que ahora usa el codigo fuente
  # (ver R/greenreg_modelo_ar.R): fitdf = p; lag = max(fitdf+1, min(10, round(n/5))), <= n-1.
  n_res    <- length(residuals(ar_ref))
  fitdf_lb <- 2
  lag_lb   <- max(fitdf_lb + 1, min(10, round(n_res / 5)))
  lag_lb   <- min(lag_lb, n_res - 1)
  lb_fitdf <- Box.test(residuals(ar_ref), lag = lag_lb, type = "Ljung-Box", fitdf = fitdf_lb)$p.value

  expect_true(grepl(sprintf("%.4f", lb_fitdf), sa, fixed = TRUE))
  # Y declara el lag/fitdf usados, para que el resultado sea auditable.
  expect_true(grepl("lag|rezago", sa, ignore.case = TRUE))
  expect_true(grepl(sprintf("fitdf=%d", fitdf_lb), sa, fixed = TRUE))
})

test_that("greenreg_modelo_ma: el mensaje de raices se refiere a raices INVERSAS", {
  b <- greenreg_modelo_ma(x, q = 1, include_mean = TRUE)
  expect_true(grepl("inversas", salida(b), ignore.case = TRUE))
})

test_that("greenreg_modelo_arma: AIC impreso = AIC(arima)", {
  cc <- greenreg_modelo_arma(x, p = 1, q = 1, include_mean = TRUE)
  arma_ref <- arima(x, order = c(1, 0, 1))
  expect_true(grepl(sprintf("%.2f", AIC(arma_ref)), salida(cc), fixed = TRUE))
})

test_that("greenreg_modelo_arima: reporta drift/deriva y media de residuos cercana a 0 (C-10)", {
  d <- greenreg_modelo_arima(y, p = 1, d = 1, q = 1, include_drift = TRUE); sd_ <- salida(d)
  expect_true(grepl("drift|deriva|Deriva", sd_, ignore.case = TRUE))
  d_sin <- greenreg_modelo_arima(y, p = 1, d = 1, q = 1); sd_sin <- salida(d_sin)
  expect_true(abs(mean(residuals(d_sin$modelo))) < 0.02 || grepl("ALERTA|ADVERTENCIA", sd_sin))
})

test_that("greenreg_modelo_ar: acepta objeto ts ademas de vector", {
  expect_s3_class(greenreg_modelo_ar(ts(x, frequency = 12), p = 1), "ar_greenreg")
})

test_that("greenreg_modelo_ar: Shapiro impreso = shapiro.test(residuos)", {
  a <- greenreg_modelo_ar(x, p = 2); ar_ref <- arima(x, order = c(2, 0, 0))
  expect_true(grepl(sprintf("%.4f", shapiro.test(residuals(ar_ref))$p.value), salida(a), fixed = TRUE))
})

test_that("C-05: greenreg_modelo_ar no aborta con n pequenio (Ljung-Box no calculable en vez de error)", {
  set.seed(2026)
  expect_no_error(print(greenreg_modelo_ar(rnorm(5), p = 2)))
})
