#' Modelo Autorregresivo de Media Móvil (ARMA)
#'
#' Esta función ajusta un modelo hibrido ARMA(p, q) que combina la dependencia
#' de los valores pasados (AutoRegresivo) con la influencia de los errores
#' pasados (Media Móvil). Es la herramienta estándar para modelar series de
#' tiempo estacionarias con estructuras de dependencia complejas.
#'
#' @details
#' El flujo interno de la función realiza una auditoría estructural en cinco pasos:
#' \enumerate{
#'   \item **Estimación Conjunta:** Utiliza Maxima Verosimilitud (MLE) para encontrar
#'         simultaneamente los coeficientes \eqn{\phi} (AR) y \eqn{\theta} (MA).
#'   \item **Estabilidad (Parte AR):** Verifica que el componente autorregresivo no sea
#'         explosivo. Las raíces inversas del polinomio AR deben estar dentro del círculo.
#'   \item **Invertibilidad (Parte MA):** Verifica que el componente de media móvil sea
#'         univoco. Las raíces inversas del polinomio MA deben estar dentro del círculo.
#'   \item **Diagnóstico de Residuos:** Evalua si el error resultante es Ruido Blanco
#'         mediante el test de Ljung-Box, asegurando que no quede señal sin capturar.
#'   \item **Criterio de Información:** Proporciona el AIC para facilitar la comparación
#'         entre diferentes combinaciones de (p, q), buscando el equilibrio entre
#'         complejidad y precisión.
#' }
#'
#' @section Guia de las 9 Gráficas de Diagnóstico (`plot`):
#' Al ejecutar \code{plot(modelo)}, se despliega una secuencia de validación mixta:
#' \describe{
#'   \item{\strong{1. Inspección de Serie}}{Muestra la serie bruta con su media. Fundamental para confirmar visualmente la estacionariedad.}
#'   \item{\strong{2. Ajuste Mixto}}{Superpone el modelo (verde) sobre la realidad. Permite evaluar como interactuan la inercia y los shocks.}
#'   \item{\strong{3. Residuos en el Tiempo}}{Evaluación de la aleatoriedad del error. Los puntos deben formar una nube sin patrones.}
#'   \item{\strong{4. Círculo Unitario Dual}}{Gráfica critica exclusiva de ARMA. Muestra simultaneamente la Estabilidad (Puntos Rojos) e Invertibilidad (Rombos Azules).}
#'   \item{\strong{5. ACF de Residuos}}{Busca autocorrelaciones remanentes. Las barras deben ser minimas para confirmar que el modelo es óptimo.}
#'   \item{\strong{6. PACF de Residuos}}{Ayuda a identificar si el orden autorregresivo (p) es suficiente.}
#'   \item{\strong{7. Normalidad (Q-Q Plot)}}{Valida que los errores sigan una campana de Gauss, requisito para la validez de los intervalos de confianza.}
#'   \item{\strong{8. Precisión (Real vs Ajustado)}}{Muestra la correlación directa. Un ARMA exitoso suele tener una alineacion estrecha sobre la diagonal.}
#'   \item{\strong{9. Estacionariedad Móvil}}{Prueba final de estabilidad. Una franja verde horizontal confirma que el modelo no se desplaza en el tiempo.}
#' }
#'
#' @param x Para \code{greenreg_modelo_arma()}: vector numérico o serie de tiempo (\code{ts}) de entrada, sin valores \code{NA}. Para \code{print()}/\code{plot()}: objeto de clase \code{"arma_greenreg"} devuelto por \code{greenreg_modelo_arma()}.
#' @param p Orden Autorregresivo (pasados de la variable).
#' @param q Orden de Media Móvil (pasados del error).
#' @param include_mean Lógico. Incluir constante/intercepto? Por defecto \code{TRUE}.
#'
#' @return Un objeto de clase \code{"arma_greenreg"} que contiene:
#' \itemize{
#'   \item \code{modelo}: Objeto \code{arima} original (orden p, 0, q).
#'   \item \code{ordenes}: Vector con los valores de p y q.
#'   \item \code{coeficientes_tabla}: Matriz con estimaciones, errores estándar y p-valores.
#'   \item \code{estatus}: Lista con el estado de estabilidad e invertibilidad.
#'   \item \code{diagnósticos}: Resultados de Ljung-Box, Shapiro-Wilk y AIC.
#'   \item \code{data_plot}: Dataframe listo para visualización con \code{ggplot2}.
#' }
#'
#' @examples
#' data("datos_anomalia_temperatura")
#' ts_data <- datos_anomalia_temperatura$anomalia_c
#' modelo <- greenreg_modelo_arma(ts_data, p = 1, q = 1)
#' modelo
#' plot(modelo)
#'
#' @importFrom stats arima coef residuals fitted Box.test pnorm qnorm shapiro.test acf pacf filter sd
#' @import ggplot2
#' @export
greenreg_modelo_arma <- function(x, p = 1, q = 1, include_mean = TRUE) {

  # --- 1. Validacion y Ajuste ---
  if (!is.numeric(x)) stop("El argumento 'x' debe ser num\u00e9rico.")

  # Ajuste ARIMA(p, 0, q)
  modelo_base <- stats::arima(x, order = c(p, 0, q), include.mean = include_mean)

  # --- 2. Coeficientes y Estadisticas ---
  coefs <- modelo_base$coef
  se <- sqrt(diag(modelo_base$var.coef))
  z_val <- coefs / se
  p_val <- 2 * (1 - stats::pnorm(abs(z_val)))

  tabla_coefs <- cbind("Estimaci\u00f3n" = coefs, "Error Est\u00e1ndar" = se, `valor z` = z_val, `Pr(>|z|)` = p_val)

  # --- 3. Analisis de Raices (Doble Verificacion) ---

  # A) Estabilidad (Parte AR) - Raices en ROJO
  ar_coefs <- coefs[grepl("^ar", names(coefs))]
  raices_ar_inv <- NULL
  es_estable <- TRUE

  if (length(ar_coefs) > 0) {
    poly_roots_ar <- polyroot(c(1, -ar_coefs))
    raices_ar_inv <- 1 / poly_roots_ar
    es_estable <- all(Mod(raices_ar_inv) < 1)
  }

  # B) Invertibilidad (Parte MA) - Raices en AZUL
  ma_coefs <- coefs[grepl("^ma", names(coefs))]
  raices_ma_inv <- NULL
  es_invertible <- TRUE

  if (length(ma_coefs) > 0) {
    poly_roots_ma <- polyroot(c(1, ma_coefs))
    raices_ma_inv <- 1 / poly_roots_ma
    es_invertible <- all(Mod(raices_ma_inv) < 1)
  }

  # --- 4. Diagnosticos de Residuos ---
  residuos <- stats::residuals(modelo_base)

  # Ljung-Box (Ruido Blanco)
  # fitdf = parametros ARMA estimados (p+q); sin esto el p-valor queda sesgado al alza.
  # Lag documentado: min(10, n/5) (criterio de forecast::checkresiduals), acotado
  # para siempre quedar entre fitdf y n (Box.test exige lag > fitdf).
  fitdf_lb <- p + q
  n_res    <- length(residuos)
  lag_lb   <- max(fitdf_lb + 1, min(10, round(n_res / 5)))
  lag_lb   <- min(lag_lb, n_res - 1)
  lb_test  <- tryCatch(
    stats::Box.test(residuos, type = "Ljung-Box", lag = lag_lb, fitdf = fitdf_lb),
    error = function(e) list(p.value = NA)
  )

  # Shapiro-Wilk (Normalidad)
  if (length(residuos) >= 3 && length(residuos) <= 5000) {
    shapiro_test <- tryCatch(stats::shapiro.test(residuos), error = function(e) list(p.value = NA))
  } else {
    shapiro_test <- list(p.value = NA)
  }

  # --- 5. Datos para Graficas ---
  ajustados <- x - residuos
  sigma <- sqrt(modelo_base$sigma2)
  # Intervalos 95%
  lower_ci <- ajustados - 1.96 * sigma
  upper_ci <- ajustados + 1.96 * sigma

  df_plot <- data.frame(
    Tiempo = 1:length(x),
    Observado = as.numeric(x),
    Ajustado = as.numeric(ajustados),
    Lower = as.numeric(lower_ci),
    Upper = as.numeric(upper_ci),
    Residuos = as.numeric(residuos)
  )

  # Datos combinados para circulo unitario
  df_raices <- data.frame(Real = numeric(), Imaginario = numeric(), Tipo = character())
  if (!is.null(raices_ar_inv)) {
    df_raices <- rbind(df_raices, data.frame(Real = Re(raices_ar_inv), Imaginario = Im(raices_ar_inv), Tipo = "AR (Estabilidad)"))
  }
  if (!is.null(raices_ma_inv)) {
    df_raices <- rbind(df_raices, data.frame(Real = Re(raices_ma_inv), Imaginario = Im(raices_ma_inv), Tipo = "MA (Invertibilidad)"))
  }

  # --- 6. Objeto Final ---
  resultado <- list(
    modelo = modelo_base,
    ordenes = c(p = p, q = q),
    coeficientes_tabla = tabla_coefs,
    estatus = list(estable = es_estable, invertible = es_invertible),
    diagnosticos = list(lb = lb_test, lb_lag = lag_lb, lb_fitdf = fitdf_lb,
                        shapiro = shapiro_test, aic = modelo_base$aic),
    data_plot = df_plot,
    data_raices = df_raices
  )

  class(resultado) <- "arma_greenreg"
  return(resultado)
}

#' @rdname greenreg_modelo_arma
#' @param x Para \code{greenreg_modelo_arma()}: vector numérico o serie de tiempo (\code{ts}) de entrada, sin valores \code{NA}. Para \code{print()}/\code{plot()}: objeto de clase \code{"arma_greenreg"} devuelto por \code{greenreg_modelo_arma()}.
#' @param ... Argumentos adicionales.
#' @export
print.arma_greenreg <- function(x, ...) {

  # Uso de sep="" para que los numeros p y q queden pegados a las comas
  cat("\n==========================================================\n")
  cat(" Modelo Mixto ARMA(", x$ordenes['p'], ",", x$ordenes['q'], ") \n", sep="")
  cat("==========================================================\n\n")

  cat("--- 1. COEFICIENTES ESTIMADOS ---\n")
  stats::printCoefmat(as.matrix(x$coeficientes_tabla), digits = 4, signif.stars = TRUE)

  # Proteccion contra NA en AIC
  aic_val <- x$diagnosticos$aic
  aic_txt <- if(is.null(aic_val) || is.na(aic_val)) "NA" else round(aic_val, 2)
  cat("\nAIC:", aic_txt, "(Menor es mejor)\n\n")

  cat("--- 2. VERIFICACI\u00d3N ESTRUCTURAL (Ra\u00edces Unitarias) ---\n")

  # Check AR (Protegido contra NA)
  if (x$ordenes['p'] > 0) {
    est <- x$estatus$estable
    if (!is.null(est) && !is.na(est) && est) {
      cat(" Parte AR: ESTABLE (Ra\u00edces dentro del c\u00edrculo).\n")
    } else {
      cat(" Parte AR: INESTABLE (Explosivo o Indeterminado).\n")
    }
  } else {
    cat(" Parte AR: No aplica (p=0).\n")
  }

  # Check MA (Protegido contra NA)
  if (x$ordenes['q'] > 0) {
    inv <- x$estatus$invertible
    if (!is.null(inv) && !is.na(inv) && inv) {
      cat(" Parte MA: INVERTIBLE (Modelo \u00fanico).\n")
    } else {
      cat(" Parte MA: NO INVERTIBLE (Problemas de identificaci\u00f3n).\n")
    }
  } else {
    cat(" Parte MA: No aplica (q=0).\n")
  }
  cat("\n")

  cat("--- 3. DIAGN\u00d3STICO DE RESIDUOS ---\n")

  # Proteccion contra NA en Ljung-Box y sprintf
  lb_p <- x$diagnosticos$lb$p.value
  lb_txt <- if(is.null(lb_p) || is.na(lb_p)) "NA" else sprintf("%.4f", lb_p)

  cat(sprintf(" Independencia (Ljung-Box): p-valor = %s  [lag=%d, fitdf=%d]\n",
              lb_txt, x$diagnosticos$lb_lag, x$diagnosticos$lb_fitdf))

  # Evaluacion logica blindada
  if (!is.null(lb_p) && !is.na(lb_p) && lb_p > 0.05) {
    cat(" Ruido Blanco: No se rechaza que el modelo haya capturado la din\u00e1mica temporal.\n")
  } else {
    cat(" Autocorrelaci\u00f3n: El modelo ARMA actual es insuficiente o indeterminado.\n")
  }

  cat("\n Usa plot(modelo) para gr\u00e1ficas de diagn\u00f3stico del modelo ARMA.\n")
  cat("==========================================================\n")

  # Retorno obligatorio S3
  return(invisible(x))
}

#' @rdname greenreg_modelo_arma
#' @description Genera una secuencia de 9 gráficas de diagnóstico para validar el modelo ARMA.
#' @param pausar Lógico. Si \code{TRUE}, pausa con \code{readline()} entre cada
#'   gráfica (por defecto, solo en sesión interactiva).
#' @param x Para \code{greenreg_modelo_arma()}: vector numérico o serie de tiempo (\code{ts}) de entrada, sin valores \code{NA}. Para \code{print()}/\code{plot()}: objeto de clase \code{"arma_greenreg"} devuelto por \code{greenreg_modelo_arma()}.
#' @param ... Argumentos adicionales.
#' @export
plot.arma_greenreg <- function(x, ..., pausar = interactive()) {

  df <- x$data_plot
  p  <- x$ordenes['p'] # Extraemos p del objeto x
  q  <- x$ordenes['q'] # Extraemos q del objeto x

  # TEMA ESTANDAR
  mi_tema <- ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, color = "#2C3E50"),
      plot.subtitle = ggplot2::element_text(size = 11, color = "#7F8C8D"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 9, face = "italic",
                                           color = "#555555", margin = ggplot2::margin(t = 15),
                                           lineheight = 1.1),
      plot.margin = ggplot2::margin(10, 10, 20, 10),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "top"
    )

  # Usar message
  if (pausar) {
    message("Generando gr\u00e1ficas ARMA... (Presiona [Enter] para avanzar)")
  }

  # --- GRAFICA 1: SERIE ORIGINAL (VALORES VS TIEMPO) ---

  media_global <- mean(df$Observado, na.rm = TRUE)

  g1 <- ggplot2::ggplot(df, ggplot2::aes(x = Tiempo, y = Observado)) +
    ggplot2::geom_line(color = "#34495E", linewidth = 0.7) +
    ggplot2::geom_point(color = "#34495E", alpha = 0.3, size = 1) +
    ggplot2::geom_hline(yintercept = media_global, linetype = "dashed",
                        color = "#16A085", linewidth = 0.8) +
    ggplot2::labs(
      title = "1. Comportamiento de la Serie Original",
      subtitle = "Visualizaci\u00f3n de la serie temporal antes del ajuste mixto",
      x = "Tiempo (Periodos)",
      y = "Valor Observado",
      caption = paste0("C\u00d3MO LEER ESTA GR\u00c1FICA:\n",
                       " L\u00cdNEA VERDE: Promedio hist\u00f3rico de la serie.\n",
                       " OBJETIVO: El modelo ARMA asume que la serie fluctua alrededor\n",
                       "   de este nivel constante (Estacionariedad en media).\n",
                       " NOTA: Si ves una 'escalera' o tendencia, el modelo ARMA(p,q)\n",
                       "   podr\u00eda requerir una diferencia previa (Modelo ARIMA).")
    ) +
    mi_tema

  print(g1)
  if (pausar) {readline(prompt = "Gr\u00e1fica 1: Inspecci\u00f3n de los datos brutos > ")}


  # --- GRAFICA 2: AJUSTE DEL MODELO ARMA (Real vs Ajustado) ---

  g2 <- ggplot2::ggplot(df, ggplot2::aes(x = Tiempo)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower, ymax = Upper),
                         fill = "#16A085", alpha = 0.15) +
    ggplot2::geom_line(ggplot2::aes(y = Observado, color = "Observado"), linewidth = 0.7) +
    ggplot2::geom_line(ggplot2::aes(y = Ajustado, color = "Modelo ARMA"), linewidth = 0.8) +
    ggplot2::scale_color_manual(values = c("Observado" = "#34495E", "Modelo ARMA" = "#16A085")) +
    ggplot2::labs(
      title = paste0("2. Ajuste del Modelo ARMA(", p, ",", q, ")"),
      subtitle = "Comparaci\u00f3n de valores reales vs. estimaci\u00f3n mixta",
      x = "Tiempo (Periodos)",
      y = "Valor",
      color = "Referencia:",
      caption = paste0("C\u00d3MO LEER ESTA GR\u00c1FICA:\n",
                       " BIEN: La l\u00ednea verde debe 'abrazar' los movimientos de la oscura.\n",
                       " SOMBRA: Zona de confianza. El 95% de los datos deben caer aqui.\n",
                       " NOTA: El ARMA es m\u00e1s potente que el AR o MA por separado,\n",
                       "   logrando un equilibrio entre tendencia y choques aleatorios.")
    ) +
    mi_tema

  print(g2)
  if (pausar) {readline(prompt = "Gr\u00e1fica 2: Evaluaci\u00f3n visual del ajuste mixto > ")}


  # --- GRAFICA 3: RESIDUOS (ERRORES) VS TIEMPO ---

  # Stats::sd
  sd_res_arma <- stats::sd(df$Residuos, na.rm = TRUE)

  g3 <- ggplot2::ggplot(df, ggplot2::aes(x = Tiempo, y = Residuos)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = -2*sd_res_arma, ymax = 2*sd_res_arma),
                         fill = "#16A085", alpha = 0.1) +
    ggplot2::geom_line(color = "#2C3E50", linewidth = 0.5) +
    ggplot2::geom_point(color = "#16A085", alpha = 0.4, size = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = c(-2*sd_res_arma, 2*sd_res_arma),
                        linetype = "dotted", color = "#E74C3C", alpha = 0.7) +
    ggplot2::labs(
      title = "3. An\u00e1lisis de Residuos (E vs t)",
      subtitle = "Verificaci\u00f3n de Aleatoriedad y Media Cero",
      x = "Tiempo (Periodos)",
      y = "Residuo (Error)",
      caption = paste0("C\u00d3MO LEER ESTA GR\u00c1FICA:\n",
                       " BIEN: Los errores deben ser una 'nube' aleatoria centrada en cero.\n",
                       " ROJO: Puntos fuera de estas l\u00edneas son shocks no explicados (Outliers).\n",
                       " MAL: Si ves ciclos, 'olas' o tendencias, el modelo no es suficiente.\n",
                       " NOTA: El ARMA busca que el error de hoy no dependa del error de ayer.")
    ) +
    mi_tema

  print(g3)
  if (pausar) {readline(prompt = "Gr\u00e1fica 3: Inspecci\u00f3n de la calidad del error residual > ")}


  # --- GRAFICA 4: ESTRUCTURA (CIRCULO UNITARIO AR Y MA) ---

  # Defensiva: validar is.null
  if (!is.null(x$data_raices) && nrow(x$data_raices) > 0) {
    theta_circ <- seq(0, 2*pi, length.out = 200)
    df_circle_arma <- data.frame(x = cos(theta_circ), y = sin(theta_circ))

    g4 <- ggplot2::ggplot() +
      ggplot2::geom_path(data = df_circle_arma, ggplot2::aes(x, y),
                         color = "#7F8C8D", linetype = "dashed", linewidth = 0.8) +
      ggplot2::geom_polygon(data = df_circle_arma, ggplot2::aes(x, y),
                            fill = "#BDC3C7", alpha = 0.1) +
      ggplot2::geom_vline(xintercept = 0, color = "gray80") +
      ggplot2::geom_hline(yintercept = 0, color = "gray80") +
      ggplot2::geom_point(data = x$data_raices,
                          ggplot2::aes(x = Real, y = Imaginario, color = Tipo, shape = Tipo),
                          size = 4, stroke = 1.2) +
      ggplot2::scale_color_manual(values = c("AR (Estabilidad)" = "#E74C3C",
                                             "MA (Invertibilidad)" = "#2980B9")) +
      ggplot2::scale_shape_manual(values = c("AR (Estabilidad)" = 16,
                                             "MA (Invertibilidad)" = 18)) +
      ggplot2::coord_fixed(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
      ggplot2::labs(
        title = "4. Estructura de Ra\u00edces (C\u00edrculo Unitario)",
        subtitle = paste0("Diagn\u00f3stico Dual: ARMA(", p, ",", q, ")"),
        x = "Parte Real",
        y = "Parte Imaginaria",
        color = "Componente:", shape = "Componente:",
        caption = paste0("C\u00d3MO LEER ESTA GR\u00c1FICA:\n",
                         " PUNTOS ROJOS (AR): Si est\u00e1n dentro, el modelo es ESTABLE.\n",
                         " ROMBOS AZULES (MA): Si est\u00e1n dentro, el modelo es INVERTIBLE.\n",
                         " CR\u00cdTICO: Si alg\u00fan punto toca o sale del c\u00edrculo gris,\n",
                         "   los resultados estadisticos NO son confiables.")
      ) +
      mi_tema + ggplot2::theme(legend.position = "right")

    print(g4)

  } else {
    # Usar message
    message("El modelo ARMA(0,0) no tiene ra\u00edces; es ruido blanco.")
  }
  if (pausar) {readline(prompt = "Gr\u00e1fica 4: Verificaci\u00f3n de estabilidad e invertibilidad mixta > ")}


  # --- GRAFICA 5: AUTOCORRELACION DE RESIDUOS (ACF) ---

  acf_res_arma <- stats::acf(df$Residuos, plot = FALSE, lag.max = 20)
  df_acf_arma <- data.frame(
    Lag = as.numeric(acf_res_arma$lag[-1]),
    ACF = as.numeric(acf_res_arma$acf[-1])
  )

  limite_arma <- 1.96 / sqrt(nrow(df))

  df_acf_arma$Etiqueta <- ifelse(abs(df_acf_arma$ACF) > limite_arma, as.character(df_acf_arma$Lag), "")
  df_acf_arma$vjust_pos <- ifelse(df_acf_arma$ACF > 0, -0.5, 1.5)

  g5 <- ggplot2::ggplot(df_acf_arma, ggplot2::aes(x = Lag, y = ACF)) +
    ggplot2::geom_col(fill = "#34495E", width = 0.4) +
    ggplot2::geom_hline(yintercept = c(limite_arma, -limite_arma),
                        linetype = "dashed", color = "#16A085", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, color = "black") +
    ggplot2::geom_text(ggplot2::aes(label = Etiqueta, vjust = vjust_pos),
                       color = "#E74C3C", fontface = "bold", size = 3.5) +
    ggplot2::scale_x_continuous(breaks = 1:20) +
    ggplot2::labs(
      title = "5. Autocorrelaci\u00f3n de Residuos (ACF)",
      subtitle = "Quedo estructura temporal en el error?",
      x = "Lag (Retraso)",
      y = "Coeficiente ACF",
      caption = paste0("C\u00d3MO LEER ESTA GR\u00c1FICA:\n",
                       " BIEN: Si todas las barras est\u00e1n 'atrapadas' en las l\u00edneas verdes.\n",
                       " N\u00daMEROS ROJOS: Indican que lags NO son ruido blanco.\n",
                       " NOTA: El ARMA busca 'limpiar' la ACF. Si el Lag 1 o 2 es rojo,\n",
                       "   el modelo mixto a\u00fan no es \u00f3ptimo.")
    ) +
    mi_tema

  print(g5)
  if (pausar) {readline(prompt = "Gr\u00e1fica 5: Diagn\u00f3stico de independencia de los errores mixtos > ")}


  # --- GRAFICA 6: AUTOCORRELACION PARCIAL (PACF) DE RESIDUOS ---

  pacf_res_arma <- stats::pacf(df$Residuos, plot = FALSE, lag.max = 20)

  df_pacf_arma <- data.frame(
    Lag = as.numeric(pacf_res_arma$lag),
    PACF = as.numeric(pacf_res_arma$acf)
  )

  limite_arma <- 1.96 / sqrt(nrow(df))

  df_pacf_arma$Etiqueta <- ifelse(abs(df_pacf_arma$PACF) > limite_arma, as.character(df_pacf_arma$Lag), "")
  df_pacf_arma$vjust_pos <- ifelse(df_pacf_arma$PACF > 0, -0.5, 1.5)

  g6 <- ggplot2::ggplot(df_pacf_arma, ggplot2::aes(x = Lag, y = PACF)) +
    ggplot2::geom_col(fill = "#16A085", width = 0.4) +
    ggplot2::geom_hline(yintercept = c(limite_arma, -limite_arma),
                        linetype = "dashed", color = "#E67E22", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, color = "black") +
    ggplot2::geom_text(ggplot2::aes(label = Etiqueta, vjust = vjust_pos),
                       color = "#E74C3C", fontface = "bold", size = 3.5) +
    ggplot2::scale_x_continuous(breaks = 1:20) +
    ggplot2::labs(
      title = "6. Autocorrelaci\u00f3n Parcial (PACF)",
      subtitle = "Identificaci\u00f3n de dependencias directas en errores mixtos",
      x = "Lag (Retraso)",
      y = "Coeficiente PACF",
      caption = paste0("C\u00d3MO LEER ESTA GR\u00c1FICA:\n",
                       " BIEN: Ninguna barra sobresale de las l\u00edneas naranjas.\n",
                       " N\u00daMEROS ROJOS: Indican que ese lag espec\u00edfico a\u00fan tiene 'se\u00f1al'.\n",
                       " NOTA: Si el Lag 1 es significativo, el componente AR(", p, ") \n",
                       "   no est\u00e1 capturando bien la inercia inmediata.")
    ) +
    mi_tema

  print(g6)
  if (pausar) {readline(prompt = "Gr\u00e1fica 6: Auditor\u00eda final de la estructura del error > ")}


  # --- GRAFICA 7: NORMALIDAD DE RESIDUOS (Q-Q PLOT) ---

  g7 <- ggplot2::ggplot(df, ggplot2::aes(sample = Residuos)) +
    ggplot2::stat_qq(color = "#16A085", alpha = 0.6, size = 2) +
    ggplot2::stat_qq_line(color = "#2C3E50", linetype = "dashed", linewidth = 1) +
    ggplot2::labs(
      title = "7. Normalidad de los Errores (Q-Q Plot)",
      subtitle = "Comparaci\u00f3n de Cuantiles de Residuos Mixtos vs. Normal Teorica",
      x = "Cuantiles Te\u00f3ricos (Normal)",
      y = "Cuantiles Observados (Residuos ARMA)",
      caption = paste0("C\u00d3MO LEER ESTA GR\u00c1FICA:\n",
                       " BIEN: Los puntos deben seguir la l\u00ednea diagonal punteada.\n",
                       " COLAS: Si los extremos se curvan, hay eventos extremos (outliers).\n",
                       " NOTA: La normalidad asegura que el error t\u00edpico sea consistente.\n",
                       " CLAVE: Si los puntos se alejan mucho, los p-valores del modelo son dudosos.")
    ) +
    mi_tema

  print(g7)
  if (pausar) {readline(prompt = "Gr\u00e1fica 7: Verificaci\u00f3n de la distribuci\u00f3n gaussiana del error > ")}


  # --- GRAFICA 8: PRECISION DEL MODELO (REAL VS AJUSTADO) ---

  cor_arma <- stats::cor(df$Observado, df$Ajustado, use = "complete.obs")

  g8 <- ggplot2::ggplot(df, ggplot2::aes(x = Ajustado, y = Observado)) +
    ggplot2::geom_point(color = "#16A085", alpha = 0.5, size = 2) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E74C3C", linewidth = 1) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "#2C3E50", fill = "#BDC3C7", alpha = 0.2) +
    ggplot2::labs(
      title = "8. Precisi\u00f3n del Modelo (Predicho vs Real)",
      subtitle = paste0("Correlaci\u00f3n de Pearson: ", round(cor_arma, 3)),
      x = "Valor Ajustado (Modelo ARMA)",
      y = "Valor Observado (Realidad)",
      caption = paste0("C\u00d3MO LEER ESTA GR\u00c1FICA:\n",
                       " L\u00cdNEA ROJA: El ideal donde el modelo predice exactamente la realidad.\n",
                       " BIEN: Los puntos verdes deben seguir la trayectoria de la l\u00ednea roja.\n",
                       " DISPERSI\u00d3N: Si los puntos est\u00e1n muy lejos de la l\u00ednea, el modelo es pobre.\n",
                       " CLAVE: Un ARMA exitoso debe tener una correlaci\u00f3n superior a 0.7.")
    ) +
    mi_tema

  print(g8)
  if (pausar) {readline(prompt = "Gr\u00e1fica 8: An\u00e1lisis de correlaci\u00f3n y capacidad predictiva > ")}


  # --- GRAFICA 9: DIAGNOSTICO DE ESTACIONARIEDAD (MOVIL) ---

  ventana_arma <- max(5, round(nrow(df) * 0.1))

  df_est_arma <- df

  df_est_arma$Media_Movil <- as.numeric(stats::filter(df$Observado, rep(1/ventana_arma, ventana_arma), sides = 2))

  # Stats::sd
  df_est_arma$SD_Movil <- sapply(1:nrow(df), function(i) {
    idx <- max(1, i-ventana_arma):min(nrow(df), i+ventana_arma)
    stats::sd(df$Observado[idx], na.rm = TRUE)
  })

  df_final_grafica <- df_est_arma[!is.na(df_est_arma$Media_Movil), ]

  g9 <- ggplot2::ggplot(df_final_grafica, ggplot2::aes(x = Tiempo)) +
    ggplot2::geom_line(data = df, ggplot2::aes(y = Observado), color = "gray85", linewidth = 0.5) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Media_Movil - SD_Movil,
                                      ymax = Media_Movil + SD_Movil),
                         fill = "#16A085", alpha = 0.15) +
    ggplot2::geom_line(ggplot2::aes(y = Media_Movil), color = "#16A085", linewidth = 1.2) +
    ggplot2::geom_hline(yintercept = mean(df$Observado, na.rm = TRUE),
                        linetype = "dashed", color = "#E74C3C", linewidth = 0.8) +
    ggplot2::labs(
      title = "9. Diagn\u00f3stico de Estacionariedad",
      subtitle = paste0("An\u00e1lisis de Media y Varianza M\u00f3vil (Ventana k = ", ventana_arma, ")"),
      x = "Tiempo (Periodos)",
      y = "Valor de la Serie",
      caption = paste0("C\u00d3MO LEER ESTA GR\u00c1FICA:\n",
                       " L\u00cdNEA VERDE: Media local. Debe ser horizontal y estable.\n",
                       " L\u00cdNEA ROJA: Promedio global. La base del modelo.\n",
                       " BIEN: Si la franja verde mantiene un ancho constante.\n",
                       " NOTA: Se han omitido los extremos para evitar errores de calculo m\u00f3vil.")
    ) +
    mi_tema

  print(g9)

  # Usar message y devolver el objeto x
  if (pausar) {
    message("Recuerda que un buen modelo ARMA es el equilibrio entre parsimonia y precisi\u00f3n.")
    message("Diagn\u00f3stico ARMA completado.")
  }

  return(invisible(list(g1 = g1, g2 = g2, g3 = g3, g4 = g4, g5 = g5, g6 = g6, g7 = g7, g8 = g8, g9 = g9)))
}