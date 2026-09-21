#' Regresion Lineal Multiple (RLM)
#'
#' Esta funcion ajusta un modelo de regresion lineal con multiples variables predictoras
#' \eqn{y = \beta_0 + \beta_1 x_1 + \beta_2 x_2 + \dots + \beta_k x_k} mediante Minimos Cuadrados Ordinarios (OLS).
#' Incluye un robusto motor de diagnostico de Multicolinealidad (VIF) y pruebas de supuestos
#' con recomendaciones de accion para el investigador.
#'
#' @details
#' El flujo interno de la funcion se divide en seis etapas:
#' \enumerate{
#'   \item **Validacion de Datos:** Asegura que la formula sea valida y que el dataframe no contenga errores estructurales.
#'   \item **Ajuste (Engine):** Utiliza \code{stats::lm}. En RLM, el foco principal es el **R2 Ajustado**,
#'         que penaliza la inclusion de variables innecesarias que no aportan informacion real.
#'   \item **Deteccion de Multicolinealidad (VIF):** Calcula manualmente el Factor de Inflacion de
#'         la Varianza \eqn{VIF = 1 / (1 - R^2_j)} para cada predictor. Esto identifica si las variables
#'         X estan correlacionadas entre si, lo cual infla los errores estandar y vuelve inestables los coeficientes.
#'   \item **Auditoria de Supuestos:** Ejecuta de forma nativa pruebas de Normalidad (Shapiro),
#'         Homocedasticidad (Breusch-Pagan), Independencia (Durbin-Watson) y Linealidad (Rainbow Test).
#'   \item **Traduccion Automatica:** Evalua los estadisticos y genera "Notas".
#'         y sugerencias de mejora (ej. transformar variables o eliminar predictores redundantes).
#'   \item **Preparacion de Diagnostico:** Construye un dataframe optimizado para \code{ggplot2} que
#'         incluye residuos estandarizados, distancias de Cook y apalancamiento (leverage).
#' }
#'
#' @section Guia de las 10 Graficas de Diagnostico (`plot`):
#' Al ejecutar \code{plot(modelo)}, se despliega una secuencia educativa de 10 visualizaciones:
#' \describe{
#'   \item{\strong{1. Distribucion de Variables (Y, Xi)}}{Muestra histogramas y densidades de todas las variables. Permite detectar sesgos o datos atipicos antes del ajuste.}
#'   \item{\strong{2. Relacion de Variables (Y vs Xi)}}{Scatterplots individuales de la respuesta contra cada predictor con una curva LOESS para verificar la linealidad visualmente.}
#'   \item{\strong{3. Ajuste (Observados vs Ajustados)}}{Muestra la precision global. Los puntos deben seguir la diagonal de 45 grados.}
#'   \item{\strong{4. Q-Q Plot (Normalidad)}}{Valida si los residuos siguen una distribucion normal. Crucial para la validez de los P-valores.}
#'   \item{\strong{5. Homocedasticidad}}{Residuos vs. Predichos. Busca una "nube de estrellas" sin forma de embudo.}
#'   \item{\strong{6. Independencia}}{Residuos vs. Orden de recoleccion. Detecta patrones temporales o espaciales en los errores.}
#'   \item{\strong{7. Prediccion vs Realidad (Identidad)}}{Version simplificada del ajuste para evaluar la correlacion de Pearson.}
#'   \item{\strong{8. Multicolinealidad (VIF)}}{Grafica de barras exclusiva de RLM. Las barras rojas (VIF > 10) indican variables que deben ser eliminadas por redundancia.}
#'   \item{\strong{9. Distancia de Cook}}{Identifica observaciones influyentes. Barras que superan el umbral $4/n$ indican datos que deforman los resultados.}
#'   \item{\strong{10. Leverage vs Residuos}}{Cruza la "rareza" de X (Leverage) con el error en Y (Residuo). Detecta puntos que actuan como palancas sobre la recta.}
#' }
#'
#' @param formula Formula modelo.
#' @param data Dataframe con los datos.
#' @param nivel_confianza Nivel de confianza para intervalos (default 0.95).
#'
#' @return Un objeto de clase \code{"rlm_greenreg"} que contiene:
#' \itemize{
#'   \item \code{modelo}: Objeto \code{lm} original.
#'   \item \code{variables}: Lista con los nombres de Y y los predictores X.
#'   \item \code{coeficientes}: Matriz con estimaciones, p-valores e intervalos de confianza.
#'   \item \code{vif}: Vector con los valores del Factor de Inflacion de la Varianza (VIF).
#'   \item \code{supuestos}: Resultados de los tests de Shapiro-Wilk (Normalidad), Breusch-Pagan 
#'         (Homocedasticidad), Durbin-Watson (Independencia) y Rainbow (Linealidad).
#'   \item \code{notas}: Reporte interpretado con recomendaciones de accion.
#'   \item \code{data_plot}: Dataframe optimizado para visualizaciones.
#' }
#'
#' @examples
#'   data(datos_rendimiento_maiz)
#'   modelo_rlm <- rlm(rendimiento_maiz_ton_ha ~ ., data = datos_rendimiento_maiz)
#'   print(modelo_rlm)
#'   plot(modelo_rlm)
#'
#' @importFrom stats lm anova shapiro.test residuals fitted coef confint sd median cor pf model.matrix cooks.distance hatvalues rstandard
#' @importFrom lmtest bptest dwtest raintest
#' @import ggplot2
#' @export
rlm <- function(formula, data, nivel_confianza = 0.95) {

  # --- 1. Validacion y Ajuste ---
  if (!inherits(formula, "formula")) stop("El argumento 'formula' debe ser una formula (ej. y ~ x1 + x2).")
  if (!is.data.frame(data)) stop("El argumento 'data' debe ser un data.frame.")

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
  tabla_coefs <- cbind(coefs, intervalos)
  tabla_anova <- stats::anova(modelo_base)

  # --- 3. Diagnosticos de Supuestos ---
  residuos <- stats::residuals(modelo_base)

  # A) Normalidad (Shapiro)
  if (length(residuos) >= 3 && length(residuos) <= 5000) {
    p_norm <- stats::shapiro.test(residuos)$p.value
  } else {
    p_norm <- NA
  }

  # B) Homocedasticidad (Breusch-Pagan)
  if (requireNamespace("lmtest", quietly = TRUE)) {
    p_homo <- lmtest::bptest(modelo_base)$p.value
  } else {
    p_homo <- NA
  }

  # C) Independencia (Durbin-Watson)
  if (requireNamespace("lmtest", quietly = TRUE)) {
    dw_res <- lmtest::dwtest(modelo_base)
    stat_dw <- dw_res$statistic
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
    test_lineal <- lmtest::raintest(modelo_base)
    p_lineal <- test_lineal$p.value
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
  nota_norm <- if (!is.na(p_norm) && p_norm > 0.05) "OKEY: Residuos normales." else " ALERTA: Residuos NO normales."
  nota_homo <- if (!is.na(p_homo) && p_homo > 0.05) "OKEY: Varianza constante." else " ALERTA: Heterocedasticidad detectada."
  nota_indep <- if (!is.na(stat_dw) && stat_dw > 1.5 && stat_dw < 2.5) "OKEY: Independencia." else " ALERTA: Posible autocorrelacion."
  nota_lineal <- if (!is.na(p_lineal) && p_lineal > 0.05) {
    "OKEY: Relacion lineal adecuada."
  } else {
    " ALERTA: Posible No-Linealidad. SUGERENCIA: Prueba transformar variables o agregar terminos cuadraticos."
  }
  # Semaforo VIF
  msg_vif <- "No aplica (1 variable)."
  if (length(vif_vals) > 0) {
    max_vif <- max(vif_vals)
    if (max_vif < 5) {
      msg_vif <- "OKEY: No hay problemas graves de Multicolinealidad (VIF < 5)."
    } else if (max_vif < 10) {
      msg_vif <- " PRECAUCIoN: Hay cierta correlacion entre variables (5 < VIF < 10)."
    } else {
      msg_vif <- "PELIGRO: Multicolinealidad severa (VIF > 10). Tus variables se estorban entre si."
    }
  }

  # --- 5. Objeto Final ---
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
      supuestos = c(nota_norm, nota_homo, nota_indep, msg_vif)
    ),
    data_plot = data.frame(
      obs_y = data[[var_y]],
      pred = fitted(modelo_base),
      resid = residuals(modelo_base),
      index = 1:length(residuos)
    )
  )

  class(resultado) <- "rlm_greenreg"
  return(resultado)
}


#' @rdname rlm
#' @param x Objeto de clase rlm_greenreg.
#' @param ... Argumentos adicionales.
#' @export
print.rlm_greenreg <- function(x, ...) {
  
  # 0. Preparacion de la Ecuacion Dinamica
  coefs_vals <- x$coeficientes[, 1]
  nombres_x <- x$variables$x
  
  # Construccion de la cadena: Y = b0 + b1*X1 + b2*X2...
  ecuacion <- paste(round(coefs_vals[1], 4)) # Intercepto
  for(i in 1:length(nombres_x)) {
    val <- round(coefs_vals[i+1], 4)
    signo <- if(val >= 0) " + " else " - "
    ecuacion <- paste0(ecuacion, signo, abs(val), " * ", nombres_x[i])
  }
  
  cat("\n==================================================================\n")
  cat("       REGRESION LINEAL MULTIPLE (RLM)        \n")
  cat("==================================================================\n\n")
  
  # --- 1. MODELO MATEMATICO ESTIMADO ---
  cat("--- 1. MODELO MATEMATICO ESTIMADO ---\n")
  cat(" Variable Dependiente (Y): ", x$variables$y, "\n")
  cat(" Predictores (X):           ", paste(x$variables$x, collapse = ", "), "\n")
  cat(" Ecuacion Estimada:         ", x$variables$y, " =", ecuacion, "\n")
  cat(" Ajuste (R2 Ajustado):      ", round(x$resumen_estadistico$R2_adj * 100, 2), "%\n\n")
  
  # --- 2. COEFICIENTES ESTIMADOS ---
  cat("--- 2. COEFICIENTES ESTIMADOS ---\n")
  # CORRECCION: Se agrega el prefijo stats:: obligatorio para CRAN
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
  cat("--- 4. ANALISIS DE VARIANZA (ANOVA) ---\n")
  if (!is.null(x$anova)) {
    # CORRECCION: Se usa printCoefmat para blindar caracteres de significancia de la ANOVA
    stats::printCoefmat(x$anova, digits = 4, pdigits = 4, has.Pvalue = TRUE, na.print = "NA")
  } else {
    cat(" No se encontro la tabla ANOVA en el objeto.\n")
  }
  cat("\n")
  
  # --- 5. VERIFICACION DE SUPUESTOS ---
  cat("--- 5. VERIFICACION DE SUPUESTOS ---\n")
  
  limpiar_txt <- function(txt) {
    if(is.null(txt) || is.na(txt) || length(txt) == 0) return("No disponible")
    res <- strsplit(txt, ": ")[[1]]
    if(length(res) > 1) return(res[2]) else return(res[1])
  }
  nota_lin   <- limpiar_txt(x$notas$supuestos[4]) # Linealidad guardada en la posicion 4
  nota_norm  <- limpiar_txt(x$notas$supuestos[1])
  nota_homo  <- limpiar_txt(x$notas$supuestos[2])
  nota_indep <- limpiar_txt(x$notas$supuestos[3])
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
                limpiar_txt(x$notas$supuestos[length(x$notas$supuestos)])))
  }
  
  # --- 6. RESUMEN DEL MODELO ---
  cat("\n--- 6. RESUMEN DEL MODELO ---\n")
  cat("GLOBAL:    ", x$notas$global, "\n")
  cat("AJUSTE:    ", x$notas$r2, "\n")
  
  # Verificacion rapida de fallos
  hay_error <- any(grepl(" ERROR", x$notas$supuestos))
  cat("VALIDEZ:   ", if(hay_error) "PRECAUCION: Se detectaron fallas en los supuestos fundamentales."
      else "El modelo cumple con los supuestos basicos.", "\n")
  
  cat("\nUse 'plot(modelo)' para ver diagnosticos visuales.\n")
  cat("==================================================================\n\n")
  return(invisible(x))
}

#' @rdname rlm
#' @description Genera una secuencia de 10 graficas para validar el modelo.
#' @param x Objeto de clase rlm_greenreg.
#' @param ... Argumentos adicionales.
#' @importFrom stats cor coef hatvalues rstandard cooks.distance
#' @importFrom ggrepel geom_text_repel
#' @export
plot.rlm_greenreg <- function(x, ...) {
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

  if (interactive()) {
    message("\nGenerando graficas de diagnostico... (Presione [Enter] para avanzar)\n")
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
      title = "1. Distribucion de Variables (Y, Xi)",
      subtitle = "Como se distribuyen los datos originales antes del ajuste?",
      x = "Valor observado",
      y = "Densidad / Frecuencia",
      caption = paste0("INTERPRETACION:\n",
                       "BIEN: Las variables muestran distribuciones equilibradas (campanas).\n",
                       " ALERTA: Si una variable esta muy sesgada (cola larga), podria afectar los coeficientes.\n",
                       " NOTA: La linea roja ayuda a identificar la forma de la distribucion.")
    ) +
    mi_tema +
    ggplot2::theme(legend.position = "none")
  
  print(g1)
  if (interactive()) {readline(prompt = "Grafica 1/10 (Distribucion Y, Xi) > ")}

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
      title = "2. Relacion de Variables (Y frente a Predictoras)",
      subtitle = paste("Visualizando la tendencia individual de cada X respecto a", var_y_name),
      x = "Valor de la Variable Predictora (Xi)",
      y = paste("Variable Respuesta (", var_y_name, ")"),
      caption = paste0("INTERPRETACION:\n",
                       "BIEN: Se observa una tendencia clara (lineal o curva suave) en los puntos.\n",
                       "MAL: Nubes de puntos totalmente dispersas sugieren que esa variable no aporta al modelo.\n",
                       "NOTA: La linea azul ayuda a identificar si la relacion es realmente lineal.")
    ) +
    mi_tema
  
  print(g2)
  if (interactive()) {readline(prompt = "Grafica 2/10 (Relacion Y vs Xi) > ")}

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
      subtitle = paste("Correlacion de ajuste (Pearson):", round(cor_ajuste, 4)),
      x = paste("Valores Reales (", x$variables$y, ")"),
      y = "Valores Predichos por el Modelo (Y)",
      caption = paste0("INTERPRETACION:\n",
                       "BIEN: Los puntos se agrupan estrechamente sobre la linea roja diagonal.\n",
                       "MAL: Nube muy dispersa o con forma curvada. El modelo no esta capturando la realidad.\n",
                       "CLAVE: La linea roja representa la perfeccion (Predicho = Real).")
    ) +
    mi_tema

  print(g3)
  if (interactive()) {readline(prompt = "Grafica 3/10 (Observados vs Ajustados) > ")}

  # ==============================================================================
  # --- GRAFICA 4: NORMALIDAD QQ ---
  # ==============================================================================
  g4 <- ggplot2::ggplot(df, ggplot2::aes(sample = resid)) +
    ggplot2::geom_qq(color = "#3498DB") +
    ggplot2::geom_qq_line(color = "#E74C3C", linetype = "dashed") +
    ggplot2::labs(
      title = "4. Grafico Q-Q (Validacion de Normalidad)",
      subtitle = "Los puntos deben seguir la linea diagonal solida.",
      x = "Cuantiles Teoricos",
      y = "Cuantiles de los Residuos",
      caption = "INTERPRETACION:\n BIEN: Los puntos estan pegados a la linea (como hormigas en fila).\n CUIDADO: Si los puntos de los extremos se curvan o alejan mucho,\ntus pruebas de significancia (P-Values) podrian ser poco fiables."
    ) + mi_tema
  
  print(g4) # CORRECCION CRAN: Se ejecuta siempre fuera del bloque condicional
  if (interactive()) {readline(prompt = "Grafica 4/10 (NORMALIDAD) > ")}

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
      caption = "INTERPRETACION:\n BIEN: Una nube rectangular sin forma definida.\n MAL (Embudo): Si los puntos se abren como un abanico <, la varianza cambia (error grave).\n MAL (Curva): Si ves una forma de 'U', te falta elevar una variable al cuadrado."
    ) + mi_tema
  
  print(g5)
  if (interactive()) {readline(prompt = "Grafica 5/10 (HOMOCEDASTICIDAD) > ")}

  # ==============================================================================
  # --- GRAFICA 6: INDEPENDENCIA ---
  # ==============================================================================
  g6 <- ggplot2::ggplot(df, ggplot2::aes(x = index, y = resid)) +
    ggplot2::geom_line(color = "#95A5A6") +
    ggplot2::geom_point(color = "#D35400") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(
      title = "6. Independencia de los Errores",
      subtitle = "Ves algun patron repetitivo u olas en los datos?",
      caption = "INTERPRETACION:\n BIEN: Ruido puro, picos arriba y abajo sin orden aparente.\n MAL: Si ves olas suaves, zig-zags constantes o grupos, \nsignifica que el pasado afecta al futuro (falta una variable temporal)."
    ) + mi_tema
  
  print(g6)
  if (interactive()) {readline(prompt = "Grafica 6/10 (INDEPENDENCIA) > ")}

  # ==============================================================================
  # --- GRAFICA 7: OBSERVADOS VS PREDICHOS ---
  # ==============================================================================
  g7 <- ggplot2::ggplot(df, ggplot2::aes(x = obs_y, y = pred)) +
    ggplot2::geom_point(color = "#2980B9", alpha = 0.6, size = 2.5) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
    ggplot2::labs(
      title = "7. Prediccion vs Realidad (Precision)",
      subtitle = "Lo ideal es que los puntos caigan sobre la linea roja discontinua.",
      x = "Valor Real (Observado)",
      y = "Valor Predicho por el Modelo",
      caption = "INTERPRETACION:\n BIEN: Los puntos forman una linea diagonal estrecha (45 grados).\n MAL: Si ves una nube redonda dispersa, el modelo no esta prediciendo nada util."
    ) + mi_tema
  
  print(g7)
  if (interactive()) {readline(prompt = "Grafica 7/10(OBSERVADOS VS PREDICHOS)  > ")}

  # ==============================================================================
  # --- GRAFICA 8: VIF (MULTICOLINEALIDAD) ---
  # ==============================================================================
  if (length(x$vif) > 0) {
    df_vif <- data.frame(Variable = names(x$vif), VIF = x$vif)
    g8 <- ggplot2::ggplot(df_vif, ggplot2::aes(x = Variable, y = VIF, fill = VIF > 10)) +
      ggplot2::geom_col() +
      ggplot2::geom_hline(yintercept = 10, color = "red", linetype = "dashed") +
      ggplot2::scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#2ECC71"), guide = "none") +
      ggplot2::labs(
        title = "8. Multicolinealidad (Inflacion de Varianza)",
        subtitle = "Variables que dicen lo mismo (Redundancia).",
        y = "Valor VIF",
        x = "Variables Predictoras",
        caption = "INTERPRETACION:\n BIEN (Verde): La variable aporta informacion unica.\n MAL (Rojo): VIF > 10. Esta variable es redundante (es combinacion de otras).\n ACCION: Elimina las variables con barra roja para estabilizar el modelo."
      ) + mi_tema
    
    print(g8)
    if (interactive()) {readline(prompt = "Grafica 8/10 (MULTICOLINEALIDAD) > ")}
  } else {
    if (interactive()) message("Grafica de VIF no generada (solo hay 1 variable predictora).")
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
      subtitle = "Algun dato esta forzando al modelo a cambiar?",
      x = "Numero de Observacion (Fila)",
      y = "Influencia (Cook)",
      caption = "INTERPRETACION:\n ALERTA: Las barras que cruzan la linea roja son datos sospechosos.\n ACCION: Revisa esas filas en tu Excel. Son errores de dedo o casos especiales?\nSi son errores, eliminalos y corre el modelo de nuevo."
    ) + mi_tema

  print(g9)
  if (interactive()) {readline(prompt = "Grafica 9/10 (INFLUENCIA) > ")}

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
    
    # CORRECCION: Usamos ggrepel ya que se encuentra declarado en las dependencias de arriba
    ggrepel::geom_text_repel(ggplot2::aes(label = Label), size = 3, fontface = "bold") +
    
    ggplot2::labs(
      title = "10. Diagnostico de Influencia (Leverage vs Residuos)",
      subtitle = "Identificacion de puntos con combinaciones de variables extremas",
      x = "Apalancamiento (Leverage)",
      y = "Residuos Estandarizados",
      caption = paste0("INTERPRETACION:\n",
                       " EJE X: Puntos a la derecha de la linea azul tienen valores de X inusuales.\n",
                       " EJE Y: Puntos fuera de las lineas rojas (2, -2) son Outliers.\n",
                       " PELIGRO: Si un punto esta muy a la derecha Y muy lejos del 0 en Y, es una 'palanca' que deforma el modelo.")
    ) +
    mi_tema +
    ggplot2::theme(legend.position = "none")

  print(g10)
  if (interactive()) {readline(prompt = "Grafica 10/10 (LEVERAGE) > ")}

  if (interactive()) {
    message("\nGraficas de diagnostico completadas.\n")
  }
  
  # CORRECCION CRAN: Devolvemos el objeto original visible/invisible segun el estandar S3
  return(invisible(x))
}
