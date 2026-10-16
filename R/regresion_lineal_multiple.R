#' Regresión Lineal Multiple (RLM)
#'
#' Esta función ajusta un modelo de regresión lineal con multiples variables predictoras
#' \eqn{y = \beta_0 + \beta_1 x_1 + \beta_2 x_2 + \dots + \beta_k x_k} mediante Minimos Cuadrados Ordinarios (OLS).
#' Incluye un robusto motor de diagnóstico de Multicolinealidad (VIF) y pruebas de supuestos
#' con recomendaciones de acción para el investigador.
#'
#' @details
#' El flujo interno de la función se divide en seis etapas:
#' \enumerate{
#'   \item **Validación de Datos:** Asegura que la fórmula sea valida y que el dataframe no contenga errores estructurales.
#'   \item **Ajuste (Engine):** Utiliza \code{stats::lm}. En RLM, el foco principal es el **R2 Ajustado**,
#'         que penaliza la inclusion de variables innecesarias que no aportan información real.
#'   \item **Detección de Multicolinealidad (VIF):** Calcula manualmente el Factor de Inflación de
#'         la Varianza \eqn{VIF = 1 / (1 - R^2_j)} para cada predictor. Esto identifica si las variables
#'         X están correlacionadas entre si, lo cual infla los errores estándar y vuelve inestables los coeficientes.
#'   \item **Auditoría de Supuestos:** Ejecuta de forma nativa pruebas de Normalidad (Shapiro),
#'         Homocedasticidad (Breusch-Pagan), Independencia (Durbin-Watson) y Linealidad (Rainbow Test).
#'   \item **Traduccion Automática:** Evalua los estadisticos y genera "Notas".
#'         y sugerencias de mejora (ej. transformar variables o eliminar predictores redundantes).
#'   \item **Preparacion de Diagnóstico:** Construye un dataframe optimizado para \code{ggplot2} que
#'         incluye residuos estandarizados, distancias de Cook y apalancamiento (leverage).
#' }
#'
#' @section Guia de las 10 Gráficas de Diagnóstico (`plot`):
#' Al ejecutar \code{plot(modelo)}, se despliega una secuencia educativa de 10 visualizaciones:
#' \describe{
#'   \item{\strong{1. Distribución de Variables (Y, Xi)}}{Muestra histogramas y densidades de todas las variables. Permite detectar sesgos o datos atípicos antes del ajuste.}
#'   \item{\strong{2. Relación de Variables (Y vs Xi)}}{Scatterplots individuales de la respuesta contra cada predictor con una curva LOESS para verificar la linealidad visualmente.}
#'   \item{\strong{3. Ajuste (Observados vs Ajustados)}}{Muestra la precisión global. Los puntos deben seguir la diagonal de 45 grados.}
#'   \item{\strong{4. Q-Q Plot (Normalidad)}}{Valida si los residuos siguen una distribución normal. Crucial para la validez de los P-valores.}
#'   \item{\strong{5. Homocedasticidad}}{Residuos vs. Predichos. Busca una "nube de estrellas" sin forma de embudo.}
#'   \item{\strong{6. Independencia}}{Residuos vs. Orden de recolección. Detecta patrones temporales o espaciales en los errores.}
#'   \item{\strong{7. Predicción vs Realidad (Identidad)}}{Version simplificada del ajuste para evaluar la correlación de Pearson.}
#'   \item{\strong{8. Multicolinealidad (VIF)}}{Gráfica de barras exclusiva de RLM. Las barras rojas (VIF > 10) indican variables que deben ser eliminadas por redundancia.}
#'   \item{\strong{9. Distancia de Cook}}{Identifica observaciones influyentes. Barras que superan el umbral $4/n$ indican datos que deforman los resultados.}
#'   \item{\strong{10. Leverage vs Residuos}}{Cruza la "rareza" de X (Leverage) con el error en Y (Residuo). Detecta puntos que actuan como palancas sobre la recta.}
#' }
#'
#' @param formula Fórmula del modelo.
#' @param data Dataframe con los datos.
#' @param nivel_confianza Nivel de confianza para intervalos (default 0.95).
#'
#' @return Un objeto de clase \code{"rlm_greenreg"} que contiene:
#' \itemize{
#'   \item \code{modelo}: Objeto \code{lm} original.
#'   \item \code{variables}: Lista con los nombres de Y y los predictores X.
#'   \item \code{coeficientes}: Matriz con estimaciones, p-valores e intervalos de confianza.
#'   \item \code{vif}: Vector con los valores del Factor de Inflación de la Varianza (VIF).
#'   \item \code{supuestos}: Resultados de los tests de Shapiro-Wilk (Normalidad), Breusch-Pagan
#'         (Homocedasticidad), Durbin-Watson (Independencia) y Rainbow (Linealidad).
#'   \item \code{notas}: Reporte interpretado con recomendaciones de acción.
#'   \item \code{data_plot}: Dataframe optimizado para visualizaciones.
#' }
#'
#' @examples
#'   data(datos_rendimiento_maiz)
#'   modelo_rlm <- greenreg_rlm(rendimiento_maiz_ton_ha ~ ., data = datos_rendimiento_maiz)
#'   print(modelo_rlm)
#'   plot(modelo_rlm)
#'
#' @importFrom stats lm anova shapiro.test residuals fitted coef confint sd median cor pf model.matrix cooks.distance hatvalues rstandard
#' @importFrom lmtest bptest dwtest raintest
#' @import ggplot2
#' @export
greenreg_rlm <- function(formula, data, nivel_confianza = 0.95) {

  # --- 1. Validacion y Ajuste ---
  if (!inherits(formula, "formula")) stop("El argumento 'f\u00f3rmula' debe ser una f\u00f3rmula (ej. y ~ x1 + x2).")
  if (!is.data.frame(data)) stop("El argumento 'data' debe ser un data.frame.")
  if (!is.numeric(nivel_confianza) || length(nivel_confianza) != 1 ||
      nivel_confianza <= 0 || nivel_confianza >= 1) {
    stop("El argumento 'nivel_confianza' debe ser un n\u00famero entre 0 y 1 (ej. 0.95).")
  }
  if (!is.numeric(nivel_confianza) || length(nivel_confianza) != 1 ||
      nivel_confianza <= 0 || nivel_confianza >= 1) {
    stop("El argumento 'nivel_confianza' debe ser un n\u00famero entre 0 y 1 (ej. 0.95).")
  }

  modelo_base <- stats::lm(formula, data = data)
  resumen <- summary(modelo_base)

  # Identificar variables
  var_y <- names(modelo_base$model)[1]
  vars_x <- names(modelo_base$model)[-1] # Todas las predictoras

  # --- 2. Estadisticos Clave ---
  # En RLM, el R2 Ajustado es el rey
  r2_adj <- resumen$adj.r.squared
  sigma <- resumen$sigma
  f_statistic <- resumen$fstatistic

  # P-valor global
  if (!is.null(f_statistic)) {
    p_valor_modelo <- pf(f_statistic[1], f_statistic[2], f_statistic[3], lower.tail = FALSE)
  } else {
    p_valor_modelo <- NA
  }

  coefs <- resumen$coefficients
  intervalos <- stats::confint(modelo_base, level = nivel_confianza)

  # Colinealidad perfecta: lm() deja NA en los coeficientes aliased; summary()
  # y confint() pueden reportar distinto numero de filas para ese caso, por
  # eso alineamos por nombre de fila en vez de pegar (cbind) a ciegas.
  colinealidad_perfecta <- anyNA(stats::coef(modelo_base))
  if (colinealidad_perfecta) {
    vars_aliased <- names(which(is.na(stats::coef(modelo_base))))
    warning(sprintf(
      "Colinealidad perfecta detectada: no se pudieron estimar los coeficientes de: %s (quedan excluidos del modelo).",
      paste(vars_aliased, collapse = ", ")
    ))
  }
  filas_comunes <- intersect(rownames(coefs), rownames(intervalos))
  tabla_coefs <- cbind(coefs[filas_comunes, , drop = FALSE],
                        intervalos[filas_comunes, , drop = FALSE])
  tabla_anova <- stats::anova(modelo_base)

  # --- 3. Diagnosticos de Supuestos ---
  residuos <- stats::residuals(modelo_base)

  # A) Normalidad (Shapiro)
  if (length(residuos) >= 3 && length(residuos) <= 5000) {
    p_norm <- tryCatch(stats::shapiro.test(residuos)$p.value, error = function(e) NA)
  } else {
    p_norm <- NA
  }

  # B) Homocedasticidad (Breusch-Pagan)
  if (requireNamespace("lmtest", quietly = TRUE)) {
    p_homo <- tryCatch(lmtest::bptest(modelo_base)$p.value, error = function(e) NA)
  } else {
    p_homo <- NA
  }

  # C) Independencia (Durbin-Watson)
  if (requireNamespace("lmtest", quietly = TRUE)) {
    stat_dw <- tryCatch(lmtest::dwtest(modelo_base)$statistic, error = function(e) NA)
  } else {
    stat_dw <- NA
  }

  # D) MULTICOLINEALIDAD (VIF) - Manual robusto para no depender de 'car'
  # VIF_j = 1 / (1 - R2_j)
  # Calculamos el VIF solo si hay mas de 1 predictor
  vif_vals <- numeric()
  if (length(vars_x) > 1) {
    tryCatch({
      # Matriz de diseno sin intercepto
      X_mat <- stats::model.matrix(modelo_base)[, -1, drop = FALSE]
      if(ncol(X_mat) > 1) {
        for(i in 1:ncol(X_mat)) {
          # Regresion de X_i contra los otros X
          y_aux <- X_mat[, i]
          x_aux <- X_mat[, -i, drop = FALSE]
          r2_aux <- summary(lm(y_aux ~ x_aux))$r.squared
          vif_vals[colnames(X_mat)[i]] <- 1 / (1 - r2_aux)
        }
      }
    }, error = function(e) {
      warning("No se pudo calcular el VIF (posiblemente variables categoricas complejas).")
    })
  }
  # E) LINEALIDAD (Rainbow Test)
  if (requireNamespace("lmtest", quietly = TRUE)) {
    p_lineal <- tryCatch(lmtest::raintest(modelo_base)$p.value, error = function(e) NA)
  } else {
    p_lineal <- NA
  }
  # --- 4. Notas (Interpretacion) ---

  # R2 Ajustado
  nota_r2 <- paste0("El modelo explica el ", round(r2_adj * 100, 2),
                    "% de la variabilidad (R2 Ajustado). ",
                    ifelse(r2_adj > 0.7, "Ajuste BUENO.",
                           ifelse(r2_adj > 0.4, "Ajuste MODERADO.", "Ajuste BAJO.")))

  # Nota Global
  nota_global <- if (!is.na(p_valor_modelo) && p_valor_modelo < 0.05) {
    "BIEN: El modelo en conjunto es SIGNIFICATIVO (al menos una variable ayuda a predecir)."
  } else {
    "MAL: El modelo en conjunto NO es significativo (ninguna variable parece servir)."
  }

  # Semaforos
  nota_norm <- if (is.na(p_norm)) {
    "No calculable."
  } else if (p_norm > 0.05) {
    "OKEY: No se rechaza la normalidad de los residuos."
  } else {
    " ALERTA: Se rechaza la normalidad de los residuos."
  }
  nota_homo <- if (is.na(p_homo)) {
    "No calculable (n insuficiente para el test)."
  } else if (p_homo > 0.05) {
    "OKEY: No se rechaza la homocedasticidad."
  } else {
    " ALERTA: Se rechaza la homocedasticidad; hay heterocedasticidad."
  }
  nota_indep <- if (is.na(stat_dw)) {
    "No calculable (n insuficiente para el test)."
  } else if (stat_dw > 1.5 && stat_dw < 2.5) {
    "OKEY: No hay evidencia de autocorrelaci\u00f3n."
  } else {
    " ALERTA: Posible autocorrelaci\u00f3n."
  }
  nota_lineal <- if (is.na(p_lineal)) {
    "No calculable (n insuficiente para el test)."
  } else if (p_lineal > 0.05) {
    "OKEY: No se rechaza la linealidad."
  } else {
    " ALERTA: Se rechaza la linealidad. SUGERENCIA: Prueba transformar variables o agregar terminos cuadraticos."
  }
  # Semaforo VIF
  msg_vif <- "No aplica (1 variable)."
  if (length(vif_vals) > 0) {
    max_vif <- max(vif_vals)
    if (max_vif < 5) {
      msg_vif <- "OKEY: No hay problemas graves de Multicolinealidad (VIF < 5)."
    } else if (max_vif < 10) {
      msg_vif <- " Precauci\u00f3n: Hay cierta correlaci\u00f3n entre variables (5 < VIF < 10)."
    } else {
      msg_vif <- "PELIGRO: Multicolinealidad severa (VIF > 10). Tus variables se estorban entre si."
    }
  }

  # --- 5. Objeto Final ---
  # Filas realmente usadas por lm() (excluye NA); fuente de obs_y para data_plot,
  # alineada con fitted()/residuals() que ya excluyen esas filas.
  mf_completo <- stats::model.frame(modelo_base)
  n_omitidas  <- nrow(data) - nrow(mf_completo)

  resultado <- list(
    modelo = modelo_base,
    variables = list(y = var_y, x = vars_x),
    coeficientes = tabla_coefs,
    anova = tabla_anova,
    vif = vif_vals,
    resumen_estadistico = list(R2_adj = r2_adj, P_valor = p_valor_modelo),
    supuestos = list(norm = p_norm, homo = p_homo, dw = stat_dw, lin = p_lineal),
    notas = list(
      r2 = nota_r2,
      global = nota_global,
      supuestos = c(norm = nota_norm, homo = nota_homo, indep = nota_indep,
                    lineal = nota_lineal, vif = msg_vif),
      na_omitidas = n_omitidas,
      colinealidad_perfecta = colinealidad_perfecta
    ),
    data_plot = data.frame(
      obs_y = mf_completo[[var_y]],
      pred = fitted(modelo_base),
      resid = residuals(modelo_base),
      index = 1:length(residuos)
    )
  )

  class(resultado) <- "rlm_greenreg"
  return(resultado)
}


#' @rdname greenreg_rlm
#' @param x Objeto de clase rlm_greenreg.
#' @param ... Argumentos adicionales.
#' @export
print.rlm_greenreg <- function(x, ...) {

  # 0. Preparacion de la Ecuacion Dinamica
  coefs_vals <- x$coeficientes[, 1]
  nombres_x <- x$variables$x

  # Construccion de la cadena: Y = b0 + b1*X1 + b2*X2...
  # Se indexa por NOMBRE (no por posicion, i+1): con colinealidad perfecta
  # (C-04) algunas variables quedan excluidas de x$coeficientes, y su
  # posicion ya no coincide con la de nombres_x.
  nombres_coef <- names(coefs_vals)
  nombres_x_presentes <- intersect(nombres_x, nombres_coef)
  intercepto <- if ("(Intercept)" %in% nombres_coef) round(coefs_vals[["(Intercept)"]], 4) else NA
  ecuacion <- paste(intercepto)
  for (nombre in nombres_x_presentes) {
    val <- round(coefs_vals[[nombre]], 4)
    signo <- if (!is.na(val) && val >= 0) " + " else " - "
    ecuacion <- paste0(ecuacion, signo, abs(val), " * ", nombre)
  }
  nombres_excluidos <- setdiff(nombres_x, nombres_x_presentes)
  if (length(nombres_excluidos) > 0) {
    ecuacion <- paste0(ecuacion, "  [excluida(s) por colinealidad perfecta: ",
                       paste(nombres_excluidos, collapse = ", "), "]")
  }

  cat("\n==================================================================\n")
  cat("       REGRESI\u00d3N LINEAL MULTIPLE (RLM)        \n")
  cat("==================================================================\n\n")

  # --- 1. MODELO MATEMATICO ESTIMADO ---
  cat("--- 1. MODELO MATEM\u00c1TICO ESTIMADO ---\n")
  cat(" Variable Dependiente (Y): ", x$variables$y, "\n")
  cat(" Predictores (X):           ", paste(x$variables$x, collapse = ", "), "\n")
  cat(" Ecuaci\u00f3n Estimada:         ", x$variables$y, " =", ecuacion, "\n")
  cat(" Ajuste (R2 Ajustado):      ", round(x$resumen_estadistico$R2_adj * 100, 2), "%\n\n")

  # --- 2. COEFICIENTES ESTIMADOS ---
  cat("--- 2. COEFICIENTES ESTIMADOS ---\n")
  # Se agrega el prefijo stats:: obligatorio para CRAN
  stats::printCoefmat(x$coeficientes[, 1:4], digits = 4, signif.stars = TRUE, has.Pvalue = TRUE)

  # Si existe VIF, lo mostramos aqui mismo para contextualizar los coeficientes
  if (length(x$vif) > 0) {
    cat("\nMulticolinealidad (VIF):\n")
    print(round(x$vif, 2))
  }
  cat("\n")

  # --- 3. INTERVALOS DE CONFIANZA ---
  cat("--- 3. INTERVALOS DE CONFIANZA ---\n")
  print(round(x$coeficientes[, 5:6], 4))
  cat("\n")

  # --- 4. ANALISIS DE VARIANZA (ANOVA) ---
  cat("--- 4. AN\u00c1LISIS DE VARIANZA (ANOVA) ---\n")
  if (!is.null(x$anova)) {
    # Se usa printCoefmat para blindar caracteres de significancia de la ANOVA
    stats::printCoefmat(x$anova, digits = 4, pdigits = 4, has.Pvalue = TRUE, na.print = "NA")
  } else {
    cat(" No se encontro la tabla ANOVA en el objeto.\n")
  }
  cat("\n")

  # --- 5. VERIFICACION DE SUPUESTOS ---
  cat("--- 5. VERIFICACI\u00d3N DE SUPUESTOS ---\n")

  limpiar_txt <- function(txt) {
    if(is.null(txt) || is.na(txt) || length(txt) == 0) return("No disponible")
    res <- strsplit(txt, ": ")[[1]]
    if(length(res) > 1) return(res[2]) else return(res[1])
  }
  nota_lin   <- limpiar_txt(x$notas$supuestos[["lineal"]])
  nota_norm  <- limpiar_txt(x$notas$supuestos[["norm"]])
  nota_homo  <- limpiar_txt(x$notas$supuestos[["homo"]])
  nota_indep <- limpiar_txt(x$notas$supuestos[["indep"]])
  p_lin_txt  <- if(is.na(x$supuestos$lin))  "NA     " else sprintf("%-7.4f", x$supuestos$lin)
  p_norm_txt <- if(is.na(x$supuestos$norm)) "NA     " else sprintf("%-7.4f", x$supuestos$norm)
  p_homo_txt <- if(is.na(x$supuestos$homo)) "NA     " else sprintf("%-7.4f", x$supuestos$homo)
  dw_txt     <- if(is.na(x$supuestos$dw))   "NA      " else sprintf("%-8.2f", x$supuestos$dw)

  cat(sprintf(" Linealidad (Rainbow Test):   [p-valor: %s] -> %s\n", p_lin_txt, nota_lin))
  cat(sprintf(" Normalidad (Shapiro-Wilk):   [p-valor: %s] -> %s\n", p_norm_txt, nota_norm))
  cat(sprintf(" Homocedasticidad (B-Pagan):  [p-valor: %s] -> %s\n", p_homo_txt, nota_homo))
  cat(sprintf(" Independencia (D-Watson):    [Estad.: %s] -> %s\n", dw_txt, nota_indep))

  # 5. VIF (Solo si hay mas de 1 X)
  if (length(x$variables$x) > 1) {
    cat(sprintf(" Multicolinealidad (VIF):      %s\n",
                limpiar_txt(x$notas$supuestos[["vif"]])))
  }

  # --- 6. RESUMEN DEL MODELO ---
  cat("\n--- 6. RESUMEN DEL MODELO ---\n")
  if (!is.null(x$notas$na_omitidas) && x$notas$na_omitidas > 0) {
    cat(sprintf(" AVISO:     Se omitieron %d fila(s) con datos faltantes (NA) para el ajuste.\n",
                x$notas$na_omitidas))
  }
  if (isTRUE(x$notas$colinealidad_perfecta)) {
    cat(" PELIGRO:   Colinealidad perfecta detectada; una o m\u00e1s variables quedaron\n")
    cat("            excluidas del modelo (revisa warnings() para ver cuales).\n")
  }
  cat("GLOBAL:    ", x$notas$global, "\n")
  cat("AJUSTE:    ", x$notas$r2, "\n")

  # Verificacion rapida de fallos
  hay_error <- any(grepl(" ERROR", x$notas$supuestos))
  cat("VALIDEZ:   ", if(hay_error) "PRECAUCI\u00d3N: Se detectaron fallas en los supuestos fundamentales."
      else "El modelo cumple con los supuestos b\u00e1sicos.", "\n")

  cat("\nUse 'plot(modelo)' para ver diagn\u00f3sticos visuales.\n")
  cat("==================================================================\n\n")
  return(invisible(x))
}

#' @rdname greenreg_rlm
#' @description Genera una secuencia de 10 gráficas para validar el modelo.
#' @param x Objeto de clase rlm_greenreg.
#' @param ... Argumentos adicionales.
#' @param pausar Lógico. Si \code{TRUE}, pausa con \code{readline()} entre cada
#'   gráfica (por defecto, solo en sesión interactiva). Con \code{FALSE} genera
#'   todas las gráficas sin pausas, útil para \code{ggsave()}, \code{patchwork} o R Markdown.
#' @importFrom stats cor coef hatvalues rstandard cooks.distance
#' @importFrom ggrepel geom_text_repel
#' @export
plot.rlm_greenreg <- function(x, ..., pausar = interactive()) {
  df <- x$data_plot
  var_y <- x$variables$y

  # Tema base
  mi_tema <- ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, color = "#2C3E50"),
      plot.subtitle = ggplot2::element_text(size = 11, color = "#34495E"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 10, face = "italic", color = "#555555", margin = ggplot2::margin(t = 10)),
      axis.title = ggplot2::element_text(face = "bold")
    )

  if (pausar) {
    message("\nGenerando gr\u00e1ficas de diagn\u00f3stico... (Presione [Enter] para avanzar)\n")
  }

  # ==============================================================================
  # --- GRAFICA 1: HISTOGRAMAS (Y, Xi) (DISTRIBUCION) ---
  # ==============================================================================
  datos_raw <- x$modelo$model
  df_dist <- data.frame(
    Valor = as.numeric(unlist(datos_raw)),
    Variable = rep(names(datos_raw), each = nrow(datos_raw))
  )

  g1 <- ggplot2::ggplot(df_dist, ggplot2::aes(x = Valor, fill = Variable)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            bins = 15, color = "white", alpha = 0.7) +
    ggplot2::geom_density(color = "#E74C3C", linewidth = 1) +
    ggplot2::facet_wrap(~Variable, scales = "free") +
    ggplot2::scale_fill_viridis_d(option = "mako", begin = 0.3, end = 0.8) +
    ggplot2::labs(
      title = "1. Distribuci\u00f3n de Variables (Y, Xi)",
      subtitle = "\u00bfC\u00f3mo se distribuyen los datos originales antes del ajuste?",
      x = "Valor observado",
      y = "Densidad / Frecuencia",
      caption = paste0("C\u00d3MO LEER ESTA GR\u00c1FICA:\n",
                       "BIEN: Las variables muestran distribuciones equilibradas (campanas).\n",
                       " ALERTA: Si una variable est\u00e1 muy sesgada (cola larga), podr\u00eda afectar los coeficientes.\n",
                       " NOTA: La l\u00ednea roja ayuda a identificar la forma de la distribuci\u00f3n.")
    ) +
    mi_tema +
    ggplot2::theme(legend.position = "none")

  print(g1)
  if (pausar) {readline(prompt = "Gr\u00e1fica 1/10 (Distribuci\u00f3n Y, Xi) > ")}

  # ==============================================================================
  # --- GRAFICA 2: RELACION DE VARIABLES (Y vs Xi) ---
  # ==============================================================================
  datos_raw <- x$modelo$model
  var_y_name <- names(datos_raw)[1]

  df_relacion <- data.frame()
  for(col in names(datos_raw)[-1]) {
    temp <- data.frame(
      Y = datos_raw[[1]],
      X_val = datos_raw[[col]],
      X_name = col
    )
    df_relacion <- rbind(df_relacion, temp)
  }

  g2 <- ggplot2::ggplot(df_relacion, ggplot2::aes(x = X_val, y = Y)) +
    ggplot2::geom_point(color = "#2C3E50", size = 2, alpha = 0.5) +
    ggplot2::geom_smooth(method = "loess", color = "#3498DB", fill = "#3498DB", alpha = 0.2, linewidth = 1) +
    ggplot2::facet_wrap(~X_name, scales = "free_x") +
    ggplot2::labs(
      title = "2. Relaci\u00f3n de Variables (Y frente a Predictoras)",
      subtitle = paste("Visualizando la tendencia individual de cada X respecto a", var_y_name),
      x = "Valor de la Variable Predictora (Xi)",
      y = paste("Variable Respuesta (", var_y_name, ")"),
      caption = paste0("C\u00d3MO LEER ESTA GR\u00c1FICA:\n",
                       "BIEN: Se observa una tendencia clara (lineal o curva suave) en los puntos.\n",
                       "MAL: Nubes de puntos totalmente dispersas sugieren que esa variable no aporta al modelo.\n",
                       "NOTA: La l\u00ednea azul ayuda a identificar si la relaci\u00f3n es realmente lineal.")
    ) +
    mi_tema

  print(g2)
  if (pausar) {readline(prompt = "Gr\u00e1fica 2/10 (Relaci\u00f3n Y vs Xi) > ")}

  # ==============================================================================
  # --- GRAFICA 3: AJUSTE (OBSERVADOS VS AJUSTADOS) ---
  # ==============================================================================
  cor_ajuste <- stats::cor(df$obs_y, df$pred)

  g3 <- ggplot2::ggplot(df, ggplot2::aes(x = obs_y, y = pred)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E74C3C", linewidth = 1) +
    ggplot2::geom_point(color = "#2C3E50", size = 3, alpha = 0.6) +
    ggplot2::geom_smooth(method = "lm", color = "#3498DB", fill = "#3498DB", alpha = 0.1, linetype = "dotted") +
    ggplot2::labs(
      title = "3. Ajuste: Valores Observados vs. Ajustados",
      subtitle = paste("Correlaci\u00f3n de ajuste (Pearson):", round(cor_ajuste, 4)),
      x = paste("Valores Reales (", x$variables$y, ")"),
      y = "Valores Predichos por el Modelo (Y)",
      caption = paste0("C\u00d3MO LEER ESTA GR\u00c1FICA:\n",
                       "BIEN: Los puntos se agrupan estrechamente sobre la l\u00ednea roja diagonal.\n",
                       "MAL: Nube muy dispersa o con forma curvada. El modelo no est\u00e1 capturando la realidad.\n",
                       "CLAVE: La l\u00ednea roja representa la perfeccion (Predicho = Real).")
    ) +
    mi_tema

  print(g3)
  if (pausar) {readline(prompt = "Gr\u00e1fica 3/10 (Observados vs Ajustados) > ")}

  # ==============================================================================
  # --- GRAFICA 4: NORMALIDAD QQ ---
  # ==============================================================================
  g4 <- ggplot2::ggplot(df, ggplot2::aes(sample = resid)) +
    ggplot2::geom_qq(color = "#3498DB") +
    ggplot2::geom_qq_line(color = "#E74C3C", linetype = "dashed") +
    ggplot2::labs(
      title = "4. Gr\u00e1fico Q-Q (Validaci\u00f3n de Normalidad)",
      subtitle = "Los puntos deben seguir la l\u00ednea diagonal solida.",
      x = "Cuantiles Te\u00f3ricos",
      y = "Cuantiles de los Residuos",
      caption = "C\u00d3MO LEER ESTA GR\u00c1FICA:\n BIEN: Los puntos est\u00e1n pegados a la l\u00ednea (como hormigas en fila).\n CUIDADO: Si los puntos de los extremos se curvan o alejan mucho,\ntus pruebas de significancia (P-Values) podr\u00edan ser poco fiables."
    ) + mi_tema

  print(g4) # Se ejecuta siempre, fuera del bloque condicional
  if (pausar) {readline(prompt = "Gr\u00e1fica 4/10 (NORMALIDAD) > ")}

  # ==============================================================================
  # --- GRAFICA 5: HOMOCEDASTICIDAD ---
  # ==============================================================================
  g5 <- ggplot2::ggplot(df, ggplot2::aes(x = pred, y = resid)) +
    ggplot2::geom_point(color = "#8E44AD", alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(
      title = "5. Homocedasticidad (Varianza Constante)",
      subtitle = "Busca una 'nube aleatoria' (cielo estrellado).",
      x = "Valores Predichos",
      y = "Errores (Residuos)",
      caption = "C\u00d3MO LEER ESTA GR\u00c1FICA:\n BIEN: Una nube rectangular sin forma definida.\n MAL (Embudo): Si los puntos se abren como un abanico <, la varianza cambia (error grave).\n MAL (Curva): Si ves una forma de 'U', te falta elevar una variable al cuadrado."
    ) + mi_tema

  print(g5)
  if (pausar) {readline(prompt = "Gr\u00e1fica 5/10 (HOMOCEDASTICIDAD) > ")}

  # ==============================================================================
  # --- GRAFICA 6: INDEPENDENCIA ---
  # ==============================================================================
  g6 <- ggplot2::ggplot(df, ggplot2::aes(x = index, y = resid)) +
    ggplot2::geom_line(color = "#95A5A6") +
    ggplot2::geom_point(color = "#D35400") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(
      title = "6. Independencia de los Errores",
      subtitle = "Ves alg\u00fan patr\u00f3n repetitivo u olas en los datos?",
      caption = "C\u00d3MO LEER ESTA GR\u00c1FICA:\n BIEN: Ruido puro, picos arriba y abajo sin orden aparente.\n MAL: Si ves olas suaves, zig-zags constantes o grupos, \nsignifica que el pasado afecta al futuro (falta una variable temporal)."
    ) + mi_tema

  print(g6)
  if (pausar) {readline(prompt = "Gr\u00e1fica 6/10 (INDEPENDENCIA) > ")}

  # ==============================================================================
  # --- GRAFICA 7: OBSERVADOS VS PREDICHOS ---
  # ==============================================================================
  g7 <- ggplot2::ggplot(df, ggplot2::aes(x = obs_y, y = pred)) +
    ggplot2::geom_point(color = "#2980B9", alpha = 0.6, size = 2.5) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
    ggplot2::labs(
      title = "7. Predicci\u00f3n vs Realidad (Precisi\u00f3n)",
      subtitle = "Lo ideal es que los puntos caigan sobre la l\u00ednea roja discontinua.",
      x = "Valor Real (Observado)",
      y = "Valor Predicho por el Modelo",
      caption = "C\u00d3MO LEER ESTA GR\u00c1FICA:\n BIEN: Los puntos forman una l\u00ednea diagonal estrecha (45 grados).\n MAL: Si ves una nube redonda dispersa, el modelo no est\u00e1 prediciendo nada \u00fatil."
    ) + mi_tema

  print(g7)
  if (pausar) {readline(prompt = "Gr\u00e1fica 7/10(OBSERVADOS VS PREDICHOS)  > ")}

  # ==============================================================================
  # --- GRAFICA 8: VIF (MULTICOLINEALIDAD) ---
  # ==============================================================================
  if (length(x$vif) > 0) {
    df_vif <- data.frame(Variable = names(x$vif), VIF = x$vif)
    vif_max <- max(x$vif)
    var_vif_max <- names(x$vif)[which.max(x$vif)]
    resultado_vif <- if (vif_max > 10) {
      sprintf(" RESULTADO: VIF m\u00e1ximo = %.1f en '%s': redundancia severa, considera eliminarla.",
              vif_max, var_vif_max)
    } else if (vif_max > 5) {
      sprintf(" RESULTADO: VIF m\u00e1ximo = %.1f en '%s': redundancia moderada, vig\u00edlala.",
              vif_max, var_vif_max)
    } else {
      sprintf(" RESULTADO: VIF m\u00e1ximo = %.1f en '%s': sin redundancia relevante.",
              vif_max, var_vif_max)
    }
    g8 <- ggplot2::ggplot(df_vif, ggplot2::aes(x = Variable, y = VIF, fill = VIF > 10)) +
      ggplot2::geom_col() +
      ggplot2::geom_hline(yintercept = 10, color = "red", linetype = "dashed") +
      ggplot2::scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#2ECC71"), guide = "none") +
      ggplot2::labs(
        title = "8. Multicolinealidad (Inflaci\u00f3n de Varianza)",
        subtitle = "Variables que dicen lo mismo (Redundancia).",
        y = "Valor VIF",
        x = "Variables Predictoras",
        caption = paste0(
          "C\u00d3MO LEER ESTA GR\u00c1FICA:\n BIEN (Verde): La variable aporta informaci\u00f3n \u00fanica.\n",
          " MAL (Rojo): VIF > 10. Esta variable es redundante (es combinaci\u00f3n de otras).\n",
          " ACCI\u00d3N: Elimina las variables con barra roja para estabilizar el modelo.\n",
          resultado_vif)
      ) + mi_tema

    print(g8)
    if (pausar) {readline(prompt = "Gr\u00e1fica 8/10 (MULTICOLINEALIDAD) > ")}
  } else {
    if (pausar) message("Gr\u00e1fica de VIF no generada (solo hay 1 variable predictora).")
  }

  # ==============================================================================
  # --- GRAFICA 9: DISTANCIA DE COOK ---
  # ==============================================================================
  cooksd <- stats::cooks.distance(x$modelo)
  n_datos <- length(cooksd)
  umbral_cook <- 4 / n_datos

  df_cook <- data.frame(Index = 1:n_datos, Cook = cooksd)
  df_cook$Label <- ifelse(df_cook$Cook > umbral_cook, as.character(df_cook$Index), "")

  g9 <- ggplot2::ggplot(df_cook, ggplot2::aes(x = Index, y = Cook)) +
    ggplot2::geom_bar(stat = "identity", fill = "#8E44AD", alpha = 0.7) +
    ggplot2::geom_hline(yintercept = umbral_cook, color = "red", linetype = "dashed") +
    ggplot2::geom_text(ggplot2::aes(label = Label), vjust = -0.5, size = 3) +
    ggplot2::labs(
      title = "9. Valores Influyentes (Distancia de Cook)",
      subtitle = "Alg\u00fan dato est\u00e1 forzando al modelo a cambiar?",
      x = "N\u00famero de Observaci\u00f3n (Fila)",
      y = "Influencia (Cook)",
      caption = "C\u00d3MO LEER ESTA GR\u00c1FICA:\n ALERTA: Las barras que cruzan la l\u00ednea roja son datos sospechosos.\n ACCI\u00d3N: Revisa esas filas en tu Excel. Son errores de dedo o casos especiales?\nSi son errores, eliminalos y corre el modelo de nuevo."
    ) + mi_tema

  print(g9)
  if (pausar) {readline(prompt = "Gr\u00e1fica 9/10 (INFLUENCIA) > ")}

  # ==============================================================================
  # --- GRAFICA 10: RESIDUOS VS APALANCAMIENTO ---
  # ==============================================================================
  lev <- stats::hatvalues(x$modelo)
  res_std <- stats::rstandard(x$modelo)
  n_datos <- length(lev)
  p_params <- length(stats::coef(x$modelo))

  umbral_lev <- 2 * (p_params / n_datos)

  df_diag <- data.frame(
    Index = 1:n_datos,
    Leverage = lev,
    ResidStd = res_std
  )

  df_diag$Label <- ifelse(df_diag$Leverage > umbral_lev | abs(df_diag$ResidStd) > 2,
                          as.character(df_diag$Index), "")

  g10 <- ggplot2::ggplot(df_diag, ggplot2::aes(x = Leverage, y = ResidStd)) +
    ggplot2::geom_point(ggplot2::aes(size = Leverage), color = "#C0392B", alpha = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "red", alpha = 0.5) +
    ggplot2::geom_vline(xintercept = umbral_lev, linetype = "dashed", color = "blue", alpha = 0.5) +

    # Usamos ggrepel ya que se encuentra declarado en las dependencias de arriba
    ggrepel::geom_text_repel(ggplot2::aes(label = Label), size = 3, fontface = "bold") +

    ggplot2::labs(
      title = "10. Diagn\u00f3stico de Influencia (Leverage vs Residuos)",
      subtitle = "Identificaci\u00f3n de puntos con combinaciones de variables extremas",
      x = "Apalancamiento (Leverage)",
      y = "Residuos Estandarizados",
      caption = paste0("C\u00d3MO LEER ESTA GR\u00c1FICA:\n",
                       " EJE X: Puntos a la derecha de la l\u00ednea azul tienen valores de X inusuales.\n",
                       " EJE Y: Puntos fuera de las l\u00edneas rojas (2, -2) son Outliers.\n",
                       " PELIGRO: Si un punto est\u00e1 muy a la derecha Y muy lejos del 0 en Y, es una 'palanca' que deforma el modelo.")
    ) +
    mi_tema +
    ggplot2::theme(legend.position = "none")

  print(g10)
  if (pausar) {readline(prompt = "Gr\u00e1fica 10/10 (LEVERAGE) > ")}

  if (pausar) {
    message("\nGraficas de diagn\u00f3stico completadas.\n")
  }

  # E-04: se devuelve la lista de graficas ggplot (en vez de invisible(x)) para
  # poder componerlas con patchwork/gridExtra, guardarlas con ggsave() o
  # reusarlas en R Markdown. g8 (VIF) solo existe con >1 variable predictora.
  graficas <- list(g1 = g1, g2 = g2, g3 = g3, g4 = g4, g5 = g5, g6 = g6, g7 = g7,
                   g9 = g9, g10 = g10)
  if (exists("g8", inherits = FALSE)) graficas$g8 <- g8
  return(invisible(graficas))
}
