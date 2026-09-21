#' Regresion Lineal Simple (RLS)
#'
#' Esta funcion ajusta un modelo de regresion lineal \eqn{y = \beta_0 + \beta_1 x} mediante
#' Minimos Cuadrados Ordinarios (OLS, por sus siglas en ingles 'Ordinary Least Squares'). 
#' A diferencia de \code{lm()} estandar, esta funcion devuelve un objeto enriquecido con 
#' interpretaciones en lenguaje natural, verificaciones automaticas de supuestos y 
#' sugerencias para el investigador.
#'
#' @details
#' El flujo interno de la funcion se divide en cinco etapas criticas:
#' \enumerate{
#'   \item **Preparacion:** Valida que los datos sean coherentes y extrae las variables de la formula.
#'   \item **Estimacion:** Calcula los coeficientes, el error estandar residual \eqn{(\sigma)} y la significancia global (F-test).
#'   \item **Inferencia:** Genera intervalos de confianza para los parametros \eqn{(\beta)} al nivel especificado.
#'   \item **Auditoria de Supuestos:** Ejecuta pruebas de Shapiro-Wilk (Normalidad), Breusch-Pagan (Homocedasticidad)
#'         y Durbin-Watson (Independencia).
#'   \item **Traduccion:** Convierte los estadisticos abstractos en "Notas" de texto que recomiendan acciones
#'         especificas (ej. aplicar logaritmos o revisar datos atipicos).
#' }
#'
#' @section Guia de las 9 Graficas de Diagnostico (\code{plot}):
#' Al ejecutar \code{plot(modelo)}, se despliega una secuencia educativa:
#' \describe{
#'   \item{\strong{1. Histogramas Exploratorios}}{Compara la distribucion de X e Y. Util para detectar sesgos iniciales.}
#'   \item{\strong{2. Relacion Original (Scatterplot)}}{Muestra la nube de puntos pura para confirmar linealidad.}
#'   \item{\strong{3. Recta de Regresion}}{Visualiza el ajuste final con bandas de confianza.}
#'   \item{\strong{4. Q-Q Plot (Normalidad)}}{Grafica los cuantiles de los residuos contra una normal teorica.}
#'   \item{\strong{5. Homocedasticidad}}{Residuos vs. Predichos para verificar varianza constante.}
#'   \item{\strong{6. Independencia}}{Residuos vs. Orden de datos para detectar autocorrelacion.}
#'   \item{\strong{7. Precision (Predicho vs Real)}}{Muestra la correlacion directa entre valores reales y estimados.}
#'   \item{\strong{8. Distancia de Cook}}{Identifica observaciones influyentes o atipicos peligrosos.}
#'   \item{\strong{9. Leverage vs Residuos}}{Diferencia entre errores en Y (Outliers) y valores extremos en X (Apalancamiento).}
#' }
#'
#' @param formula Una formula de R estandar (ej. \code{rendimiento ~ fertilizante}).
#' @param data El data.frame que contiene las variables.
#' @param nivel_confianza Nivel de confianza para los intervalos (por defecto 0.95).
#'
#' @return Un objeto de clase \code{"rls_greenreg"} que es una lista con:
#' \itemize{
#'   \item \code{modelo}: El objeto \code{lm} original.
#'   \item \code{variables}: Nombres de X e Y.
#'   \item \code{coeficientes}: Matriz con Estimaciones, Errores e Intervalos.
#'   \item \code{supuestos}: Resultados de los tests de Shapiro, BP y DW.
#'   \item \code{notas}: Diccionario de interpretaciones y sugerencias de mejora.
#'   \item \code{data_plot}: Dataframe con residuos e indices de influencia.
#' }
#'
#' @examples
#'   data(datos_rendimiento_maiz)
#'   modelo <- rls(rendimiento_maiz_ton_ha ~ precipitacion_mm, data = datos_rendimiento_maiz)
#'   print(modelo)
#'   plot(modelo)
#'
#' @importFrom stats lm anova shapiro.test residuals fitted coef confint sd median cor pf cooks.distance hatvalues rstandard
#' @importFrom lmtest bptest dwtest
#' @import ggplot2
#' @export
  rls <- function(formula, data, nivel_confianza = 0.95) {
    
    if (!inherits(formula, "formula")) stop("El argumento 'formula' debe ser una formula (ej. y ~ x).")
    if (!is.data.frame(data)) stop("El argumento 'data' debe ser un data.frame.")
    
    modelo_base <- stats::lm(formula, data = data)
    resumen <- summary(modelo_base)
    
    var_y <- names(modelo_base$model)[1]
    var_x <- names(modelo_base$model)[2]
    
    r2 <- resumen$r.squared
    r2_adj <- resumen$adj.r.squared
    sigma <- resumen$sigma
    f_statistic <- resumen$fstatistic
    
    if (!is.null(f_statistic)) {
      p_valor_modelo <- stats::pf(f_statistic[1], f_statistic[2], f_statistic[3], lower.tail = FALSE)
    } else {
      p_valor_modelo <- NA
    }
    
    coefs <- resumen$coefficients
    intervalos <- stats::confint(modelo_base, level = nivel_confianza)
    tabla_coefs <- cbind(coefs, intervalos)
    tabla_anova <- stats::anova(modelo_base)
    
    residuos <- stats::residuals(modelo_base)
    
    # A) Normalidad (Shapiro-Wilk)
    if (length(residuos) >= 3 && length(residuos) <= 5000) {
      test_shapiro <- stats::shapiro.test(residuos)
      p_norm <- test_shapiro$p.value
    } else {
      p_norm <- NA
    }
    
    # B) Homocedasticidad (Breusch-Pagan)
    p_homo <- NA
    if (requireNamespace("lmtest", quietly = TRUE)) {
      test_bp <- lmtest::bptest(modelo_base)
      p_homo <- test_bp$p.value
    }
    
    # C) Independencia (Durbin-Watson)
    stat_dw <- NA; p_indep <- NA
    if (requireNamespace("lmtest", quietly = TRUE)) {
      test_dw <- lmtest::dwtest(modelo_base)
      stat_dw <- test_dw$statistic
      p_indep <- test_dw$p.value
    }
    
    # --- Generacion de Notas ---
    nota_r2 <- paste0("El modelo explica el ", round(r2 * 100, 2),
                      "% de la variabilidad de '", var_y, "'. ",
                      ifelse(r2 > 0.7, "Es un ajuste ALTO (Bueno).",
                             ifelse(r2 > 0.4, "Es un ajuste MODERADO.", "Es un ajuste BAJO.")))
    
    beta1 <- coefs[2, 1]
    pval_beta1 <- coefs[2, 4]
    nota_beta1 <- if (pval_beta1 < 0.05) {
      paste0("La relacion es SIGNIFICATIVA. Por cada aumento de 1 unidad en '", var_x,
             "', la variable '", var_y, "' cambia en promedio ", round(beta1, 4), " unidades.")
    } else {
      paste0("La variable '", var_x, "' NO tiene un efecto estadisticamente significativo. ",
             "SUGERENCIA: Intenta recolectar mas datos o revisa si la relacion es curva.")
    }
    
    nota_norm <- if (!is.na(p_norm) && p_norm > 0.05) {
      "OKEY: Los residuos parecen normales (p > 0.05)."
    } else {
      "ALERTA: Residuos NO normales. ACCION: Prueba transformar Y o usa modelos no parametricos."
    }
    
    nota_homo <- if (!is.na(p_homo) && p_homo > 0.05) {
      "OKEY: Varianza constante (Homocedasticidad)."
    } else {
      "ALERTA: Heterocedasticidad detectada. ACCION: Se recomienda usar errores estandar robustos."
    }
    
    nota_indep <- if (!is.na(stat_dw) && stat_dw > 1.5 && stat_dw < 2.5) {
      "OKEY: Residuos independientes (DW cercano a 2)."
    } else {
      "ALERTA: Posible autocorrelacion. ACCION: Revisa si los datos son una serie temporal."
    }
    
    fallo_critico <- (!is.na(p_norm) && p_norm <= 0.05) || (!is.na(p_homo) && p_homo <= 0.05) || (!is.na(stat_dw) && (stat_dw <= 1.5 | stat_dw >= 2.5))
    nota_experto <- if(fallo_critico) "RECOMENDACION FINAL: Fallas en supuestos. Consulta a un analista." else "El modelo cumple con los supuestos basicos."
    
    resultado <- list(
      modelo = modelo_base,
      variables = list(y = var_y, x = var_x),
      info_datos = list(n = length(residuos), media_y = mean(data[[var_y]], na.rm=TRUE), sd_y = stats::sd(data[[var_y]], na.rm=TRUE)),
      anova = tabla_anova,
      coeficientes = tabla_coefs,
      resumen_estadistico = list(R2 = r2, R2_adj = r2_adj, Sigma = sigma, F_stat = f_statistic),
      supuestos = list(shapiro_p = p_norm, bp_p = p_homo, dw_stat = stat_dw),
      notas = list(r2 = nota_r2, impacto = nota_beta1, supuestos = c(nota_norm, nota_homo, nota_indep), alerta_experto = nota_experto),
      data_plot = data.frame(obs_y = data[[var_y]], obs_x = data[[var_x]], pred = stats::fitted(modelo_base), resid = residuos, index = 1:length(residuos))
    )
    
    class(resultado) <- "rls_greenreg"
    return(resultado)
  }

  #' @rdname rls
  #' @param x Objeto de clase rls_greenreg.
  #' @param ... Argumentos adicionales.
  #' @export
  print.rls_greenreg <- function(x, ...) {
    b0 <- round(x$coeficientes[1, 1], 4)
    b1 <- round(x$coeficientes[2, 1], 4)
    signo_txt <- if(b1 >= 0) " + " else " - "
    
    cat("\n==================================================================\n")
    cat("       REGRESION LINEAL SIMPLE        \n")
    cat("==================================================================\n\n")
    
    cat("\n--- 1. MODELO MATEMATICO ESTIMADO ---\n")
    cat(" Variable Dependiente (Y): ", x$variables$y, "\n")
    cat(" Variable Independiente (X):", x$variables$x, "\n")
    cat(" Ecuacion de la Recta:   ", x$variables$y, " =", b0, signo_txt, abs(b1), "*", x$variables$x, "\n")
    cat(" Bondad de Ajuste (R2):   ", round(x$resumen_estadistico$R2 * 100, 2), "%\n\n")
    
    cat("\n--- 2. COEFICIENTES ESTIMADOS ---\n")
    stats::printCoefmat(x$coeficientes[, 1:4], digits = 4, signif.stars = TRUE)
    
    # --- NUEVA SECCIÓN: TABLA ANOVA INCORPORADA ---
    cat("\n--- 3. TABLA DE ANALISIS DE VARIANZA (ANOVA) ---\n")
    if (!is.null(x$anova)) {
      # Usamos printCoefmat para controlar los digitos y evitar caracteres raros de las estrellas en ANOVA
      stats::printCoefmat(x$anova, digits = 4, pdigits = 4, has.Pvalue = TRUE, na.print = "NA")
    } else {
      cat(" No se encontro la tabla ANOVA en el objeto.\n")
    }
    
    cat("\n--- 4. VERIFICACION DE SUPUESTOS ---\n")
    
    # Correccion para p-valor de Shapiro (Evita caidas por NA)
    p_shapiro_txt <- if(is.na(x$supuestos$shapiro_p)) "NA     " else sprintf("%-7.4f", x$supuestos$shapiro_p)
    cat(sprintf(" Normalidad (Shapiro-Wilk):   [p-valor: %s] -> %s\n", p_shapiro_txt, x$notas$supuestos[1]))
    
    # Correccion para p-valor de Breusch-Pagan (Evita caidas por NA)
    p_bp_txt <- if(is.na(x$supuestos$bp_p)) "NA     " else sprintf("%-7.4f", x$supuestos$bp_p)
    cat(sprintf(" Homocedasticidad (B-Pagan):  [p-valor: %s] -> %s\n", p_bp_txt, x$notas$supuestos[2]))
    
    # Correccion para estadistico Durbin-Watson (Evita caidas por NA)
    dw_txt <- if(is.na(x$supuestos$dw_stat)) "NA      " else sprintf("%-8.2f", x$supuestos$dw_stat)
    cat(sprintf(" Independencia (D-Watson):    [Estad.: %s] -> %s\n", dw_txt, x$notas$supuestos[3]))
    
    cat("\n--- 5. RESUMEN DEL MODELO ---\n")
    cat("IMPACTO:  ", x$notas$impacto, "\n")
    cat("VALIDEZ:  ", x$notas$alerta_experto, "\n")
    
    return(invisible(x))
  }
  
  #' @rdname rls
  #' @description Genera una secuencia de 9 graficas para validar el modelo.
  #' @param x Objeto de clase rls_greenreg.
  #' @param ... Argumentos adicionales.
  #' @export
  plot.rls_greenreg <- function(x, ...) {
    
    df <- x$data_plot
    var_y <- x$variables$y
    var_x <- x$variables$x
    
    # Tema visual consistente
    mi_tema <- ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14, color = "#2C3E50"),
        plot.subtitle = ggplot2::element_text(size = 11, color = "#34495E"),
        plot.caption = ggplot2::element_text(hjust = 0, size = 10, face = "italic", color = "#555555", margin = ggplot2::margin(t = 10)),
        axis.title = ggplot2::element_text(face = "bold")
      )
    
    # CORRECCION CRAN: Usamos message() en lugar de cat()
    if (interactive()) {
      message("\nGenerando graficas de diagnostico...")
      message("Presione [Enter] en la consola para avanzar entre las 9 visualizaciones.\n")
    }
    
    # ==============================================================================
    # --- GRAFICA 1: HISTOGRAMAS EXPLORATORIOS (X, Y) ---
    # ==============================================================================
    df_long <- data.frame(
      Valor = c(df$obs_x, df$obs_y),
      Variable = c(rep(paste("Variable X:", var_x), nrow(df)),
                   rep(paste("Variable Y:", var_y), nrow(df)))
    )
    
    g1 <- ggplot2::ggplot(df_long, ggplot2::aes(x = Valor, fill = Variable)) +
      ggplot2::geom_histogram(bins = 15, color = "white", alpha = 0.8, show.legend = FALSE) +
      ggplot2::facet_wrap(~Variable, scales = "free") + 
      ggplot2::scale_fill_manual(values = c("#3498DB", "#16A085")) +
      ggplot2::labs(
        title = "1. Exploracion Inicial de Datos (Histogramas)",
        subtitle = "Como se distribuyen las variables originales antes del modelo?",
        x = "Valor de la Variable",
        y = "Frecuencia Absoluta",
        caption = paste0("INTERPRETACION:\n",
                         "BIEN: Las barras forman una figura aproximadamente simetrica (campana).\n",
                         " ALERTA: Si un histograma tiene una 'cola' muy larga hacia la derecha o izquierda,\n",
                         "   la variable esta sesgada. Esto puede afectar la precision de la recta de regresion.")
      ) +
      mi_tema
    print(g1)
    if (interactive()) {readline(prompt = "Grafica 1/9 (Histogramas X,Y) > Presione [Enter] para continuar...")}
    
    # ==============================================================================
    # --- GRAFICA 2: RELACION ORIGINAL (SCATTERPLOT) ---
    # ==============================================================================
    g2 <- ggplot2::ggplot(df, ggplot2::aes(x = obs_x, y = obs_y)) +
      ggplot2::geom_point(color = "#2C3E50", size = 3, alpha = 0.6) +
      ggplot2::labs(
        title = "2. Relacion Original (Scatterplot)",
        subtitle = paste("Explorando la correlacion visual entre", var_x, "y", var_y),
        x = var_x,
        y = var_y,
        caption = paste0("INTERPRETACION:\n",
                         "BIEN: Los puntos forman una hilera o elipse clara (ascendente o descendente).\n",
                         "MAL: Los puntos parecen una nube circular o ruidosa. Indica que no hay relacion lineal.")
      ) +
      mi_tema
    
    print(g2)
    if (interactive()) {readline(prompt = "Grafica 2/9 (Relacion Original) > ")}
    
    # ==============================================================================
    # --- GRAFICA 3: RECTA DE REGRESION (AJUSTE) ---
    # ==============================================================================
    g3 <- ggplot2::ggplot(df, ggplot2::aes(x = obs_x, y = obs_y)) +
      ggplot2::geom_point(color = "#2C3E50", size = 3, alpha = 0.4) +
      ggplot2::geom_smooth(method = "lm", color = "#E74C3C", fill = "#E74C3C", alpha = 0.2, linewidth = 1.2) +
      ggplot2::labs(
        title = "3. Ajuste del Modelo (Recta de Regresion)",
        subtitle = paste("Modelo Estimado:", var_y, "~", var_x),
        x = var_x,
        y = var_y,
        caption = paste0("INTERPRETACION:\n",
                         "BIEN: La linea roja cruza por el centro de la nube de puntos.\n",
                         "NOTA: La banda sombreada es el error de estimacion. Si es muy ancha,\n",
                         "   el modelo tiene mucha incertidumbre en esa zona.")
      ) +
      mi_tema
    
    print(g3)
    if (interactive()) {readline(prompt = "Grafica 3/9 (Ajuste del Modelo) > ")}
    
    # ==============================================================================
    # --- GRAFICA 4: NORMALIDAD (QQ-PLOT) ---
    # ==============================================================================
    g4 <- ggplot2::ggplot(df, ggplot2::aes(sample = resid)) +
      ggplot2::geom_qq(color = "#3498DB", size = 2) +
      ggplot2::geom_qq_line(color = "#E74C3C", linewidth = 1.2, linetype = "dashed") +
      ggplot2::labs(
        title = "4. Validacion de Normalidad (Q-Q Plot)",
        subtitle = "Comparacion contra una distribucion normal ideal",
        x = "Cuantiles Teoricos", y = "Cuantiles Observados",
        caption = "INTERPRETACION:\nBIEN: Los puntos azules estan pegados a la linea roja discontinua.\n MAL: Si los puntos en los extremos se curvan o alejan mucho de la linea."
      ) + mi_tema
    print(g4)
    if (interactive()) {readline(prompt = "Grafica 4/9 (Normalidad) > ")}
    
    # ==============================================================================
    # --- GRAFICA 5: HOMOCEDASTICIDAD ---
    # ==============================================================================
    g5 <- ggplot2::ggplot(df, ggplot2::aes(x = pred, y = resid)) +
      ggplot2::geom_point(color = "#8E44AD", size = 2.5, alpha = 0.7) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      ggplot2::labs(
        title = "5. Homocedasticidad (Varianza Constante)",
        subtitle = "El error es igual para valores chicos y grandes?",
        x = "Valores Predichos", y = "Residuos (Errores)",
        caption = "INTERPRETACION:\n BIEN: Nube dispersa sin forma (como cielo estrellado).\n MAL: Forma de 'Embudo' (<) o abanico. Significa que el error crece con el valor."
      ) + mi_tema
    print(g5)
    if (interactive()) {readline(prompt = "Grafica 5/9 Homocedasticidad > ")}
    
    # ==============================================================================
    # --- GRAFICA 6: INDEPENDENCIA ---
    # ==============================================================================
    g6 <- ggplot2::ggplot(df, ggplot2::aes(x = index, y = resid)) +
      ggplot2::geom_line(color = "#95A5A6", alpha = 0.5) +
      ggplot2::geom_point(color = "#D35400", size = 2) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::labs(
        title = "6. Independencia (Residuos vs Orden)",
        subtitle = "Detectando patrones en el tiempo de recoleccion",
        x = "Orden de Recoleccion", y = "Residuos",
        caption = "INTERPRETACION:\n BIEN: Ruido aleatorio (picos arriba y abajo sin orden).\n MAL: Si ves ondas suaves o tendencias, los datos no son independientes."
      ) + mi_tema
    print(g6)
    if (interactive()) { readline(prompt = "Grafica 6/9 Independencia > ")}
    
    # ==============================================================================
    # --- GRAFICA 7: PREDICCION VS REALIDAD ---
    # ==============================================================================
    correlacion <- stats::cor(df$obs_y, df$pred)
    g7 <- ggplot2::ggplot(df, ggplot2::aes(x = obs_y, y = pred)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
      ggplot2::geom_point(color = "#27AE60", size = 3, alpha = 0.7) +
      ggplot2::labs(
        title = "7. Precision: Prediccion vs Realidad",
        subtitle = paste("Correlacion:", round(correlacion, 3)),
        x = "Valor Real", y = "Valor Predicho",
        caption = "INTERPRETACION:\n BIEN: Los puntos caen sobre la linea gris punteada (Identidad).\n MAL: Nube dispersa y redonda (el modelo no predice bien)."
      ) + mi_tema
    print(g7)
    if (interactive()) {readline(prompt = "Grafica 7/9 > ")}
    
    # ==============================================================================
    # --- GRAFICA 8: DISTANCIA DE COOK ---
    # ==============================================================================
    cooksd <- stats::cooks.distance(x$modelo)
    n_datos <- length(cooksd)
    umbral_cook <- 4 / n_datos
    df_cook <- data.frame(Index = 1:n_datos, Cook = cooksd)
    df_cook$Label <- ifelse(df_cook$Cook > umbral_cook, as.character(df_cook$Index), "")
    
    g8 <- ggplot2::ggplot(df_cook, ggplot2::aes(x = Index, y = Cook)) +
      ggplot2::geom_bar(stat = "identity", fill = "#E67E22", alpha = 0.8) +
      ggplot2::geom_hline(yintercept = umbral_cook, color = "red", linetype = "dashed") +
      ggplot2::geom_text(ggplot2::aes(label = Label), vjust = -0.5, size = 3) +
      ggplot2::labs(
        title = "8. Valores Influyentes (Distancia de Cook)",
        subtitle = "Deteccion de datos atipicos peligrosos",
        x = "Numero de Observacion", y = "Influencia",
        caption = "INTERPRETACION:\n ALERTA: Las barras que cruzan la linea roja son datos sospechosos.\n ACCION: Revisa esas filas en tu Excel, podrian ser errores de captura."
      ) + mi_tema
    print(g8)
    if (interactive()) {readline(prompt = "Grafica 8/9 Influyentes > ")}
    
    # ==============================================================================
    # --- GRAFICA 9: RESIDUOS VS APALANCAMIENTO ---
    # ==============================================================================
    lev <- stats::hatvalues(x$modelo)
    res_std <- stats::rstandard(x$modelo)
    n_datos <- length(lev)
    umbral_lev <- 4 / n_datos
    
    df_diag <- data.frame(
      index = 1:n_datos,
      leverage = lev,
      resid_std = res_std
    )
    
    df_diag$Label <- ifelse(df_diag$leverage > umbral_lev | abs(df_diag$resid_std) > 2,
                            as.character(df_diag$index), "")
    
    g9 <- ggplot2::ggplot(df_diag, ggplot2::aes(x = leverage, y = resid_std)) +
      ggplot2::geom_point(ggplot2::aes(size = leverage), color = "#C0392B", alpha = 0.6) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "red", alpha = 0.5) +
      ggplot2::geom_vline(xintercept = umbral_lev, linetype = "dashed", color = "blue", alpha = 0.5)
    
    # Verificacion defensiva del paquete ggrepel
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      g9 <- g9 + ggrepel::geom_text_repel(ggplot2::aes(label = Label), size = 3.5, fontface = "bold")
    } else {
      g9 <- g9 + ggplot2::geom_text(ggplot2::aes(label = Label), size = 3.5, vjust = -0.5)
    }
    
    g9 <- g9 + 
      ggplot2::labs(
        title = "9. Diagnostico de Influencia (Leverage vs Residuos)",
        subtitle = "Identificando puntos que 'jalonean' la recta",
        x = "Apalancamiento (Leverage)",
        y = "Residuos Estandarizados",
        caption = paste0("INTERPRETACION:\n",
                         " EJE X: Puntos a la derecha de la linea azul tienen valores de X muy raros.\n",
                         " EJE Y: Puntos fuera de las lineas rojas punteadas son Outliers (errores en Y).\n",
                         " CLAVE: Si un punto esta muy a la derecha Y muy arriba/abajo, esta deformando tu modelo.")
      ) +
      mi_tema +
      ggplot2::theme(legend.position = "none")
    
    print(g9)
    
    # CORRECCION CRAN: Usamos message()
    if (interactive()) {
      message("Secuencia de graficas completada.\n")
    }
  }