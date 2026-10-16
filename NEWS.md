# GREENREG 0.1.2 (en desarrollo)

Esta version corrige los hallazgos de severidad Alta identificados en la
revision de codigo de septiembre 2026 (M.C. Israel Lerma Serna, DICIFO-UACh).

## Cambio que rompe compatibilidad (con alias obsoletos)

* Las 9 funciones exportadas se renombraron con el prefijo `greenreg_`:
  `rlm()` -> `greenreg_rlm()`, `rls()` -> `greenreg_rls()`,
  `logistico()` -> `greenreg_logistico()`, `reg_poisson()` -> `greenreg_reg_poisson()`,
  `modelo_ar()` -> `greenreg_modelo_ar()`, `modelo_ma()` -> `greenreg_modelo_ma()`,
  `modelo_arma()` -> `greenreg_modelo_arma()`, `modelo_arima()` -> `greenreg_modelo_arima()`,
  `analisis_datos()` -> `greenreg_analisis_datos()`.
  Motivo (hallazgo E-03): `rlm()` enmascaraba `MASS::rlm()` (regresion robusta)
  sin ningun aviso; el orden en que se cargaban los paquetes decidia en
  silencio que modelo se ajustaba. Los nombres antiguos se mantienen como
  alias obsoletos (`.Deprecated()`, ver `R/deprecated.R`) durante esta
  version y seguiran funcionando igual, mostrando un aviso.

## Correcciones

* `greenreg_rlm()`: la linea "Linealidad (Rainbow Test)" imprimia el mensaje del VIF
  en vez del propio; `notas$supuestos` ahora es un vector con nombres para
  evitar este tipo de desalineacion por posicion.
* `greenreg_rlm()`/`greenreg_rls()`: `data_plot` ya no aborta cuando hay `NA` en los
  predictores; ademas se avisa cuantas filas se omitieron por datos
  faltantes.
* `greenreg_logistico()`: el VIF ya no depende de `car::vif()` (se calcula
  manualmente) y ya no aborta con un solo predictor.
* `greenreg_rlm()`: la colinealidad perfecta entre predictores ya no rompe la tabla
  de coeficientes; ahora se reporta como una alerta explicita ("PELIGRO:
  Colinealidad perfecta...") en vez de detener la ejecucion.
* `greenreg_rlm()`/`greenreg_modelo_ar()`: las pruebas auxiliares (Breusch-Pagan,
  Durbin-Watson, Rainbow, Ljung-Box) ya no abortan con muestras pequenias;
  cuando no son calculables, se reporta "No calculable" en vez de un error
  interno de R.
* `greenreg_modelo_ar()`, `greenreg_modelo_ma()`, `greenreg_modelo_arma()`, `greenreg_modelo_arima()`: el test
  de Ljung-Box ahora usa `fitdf` (parametros ARMA estimados), evitando un
  p-valor sesgado al alza; el `lag` y el `fitdf` usados ahora se declaran
  explicitamente en la salida.

## Estructura del paquete

* Se eliminaron cuatro archivos duplicados en `R/` (`Modelo_*.R`) que
  contenian copias identicas, sin mantener, de `greenreg_modelo_ar.R`, `greenreg_modelo_ma.R`,
  `greenreg_modelo_arma.R` y `greenreg_modelo_arima.R`.
* Se agrego una suite de pruebas unitarias con `testthat` (`tests/testthat/`)
  basada en las comprobaciones funcionales y de estres de la revision de
  codigo.
* Se agregaron `README.md`, este `NEWS.md` y una viñeta de introduccion.

## Correcciones (severidad Media)

* `greenreg_reg_poisson()`: la ecuacion impresa usaba `round(coef, 3)`, redondeando a
  0 coeficientes pequenos pero significativos; ahora usa `signif(coef, 4)` y
  el texto dice `log(mu)` en vez de `log(Y)`.
* `greenreg_logistico()`: la tabla de Odds Ratios ya no incluye el intercepto (su OR
  no tiene interpretacion util y aparecia como `0.0000` por redondeo).
* `greenreg_logistico()`: se agregan sensibilidad, especificidad, VPP, VPN, tasa de
  no informacion y exactitud balanceada junto a la accuracy; la nota del
  Pseudo R2 (McFadden) ahora se condiciona al rango 0.2-0.4 en vez de
  afirmarlo siempre.
* `greenreg_modelo_arima()`: nuevo argumento `include_drift`; con `d >= 1` se avisa
  si la media de residuos sugiere una tendencia sin capturar.
* `greenreg_analisis_datos()`: la deteccion de variables de conteo se unifico en una
  sola funcion interna que exige enteros no negativos y una heuristica
  adicional (pocos valores distintos o media baja), evitando falsos
  positivos con variables continuas medidas en enteros (edad, cm, ingresos).
* `greenreg_analisis_datos()`: la deteccion de series de tiempo ahora acepta objetos
  `ts` ademas de columnas de fecha, y usa Ljung-Box global en vez de "algun
  rezago de la ACF supera el umbral" (que daba falsos positivos frecuentes).
* `greenreg_analisis_datos()`: la "Opcion B: modelo optimizado" ahora indica
  explicitamente si el stepwise (AIC) conservo todas las variables o cuales
  elimino, en vez de imprimir la misma formula que la Opcion A sin avisar.
* `greenreg_rlm()`/`greenreg_rls()`: se valida `nivel_confianza` (debe estar entre 0 y 1);
  `greenreg_rls()` ahora exige exactamente un predictor.

## Documentacion

* Se corrigieron dos referencias con ISBN incorrecto en `DESCRIPTION`: la
  cita de Montgomery et al. ahora corresponde a *Introduction to Linear
  Regression Analysis*, 6.a ed. (Montgomery, Peck & Vining, 2021, Wiley,
  ISBN 978-1119578727), y la de Box & Jenkins a *Time Series Analysis:
  Forecasting and Control*, 5.a ed. (Box, Jenkins, Reinsel & Ljung, 2015,
  Wiley, ISBN 978-1118675021), ambos ISBN verificados.
* `DESCRIPTION` (`Authors@R`): los cuatro miembros del comite de tesis
  listados como contribuidores (`ctb`, que implica contribucion de codigo)
  pasan a `rev` (revisor), su rol real; se confirmo su consentimiento
  explicito para aparecer en el paquete publicado.
* Se corrigieron acentos y eñes en las cadenas de texto y comentarios
  roxygen de los 9 archivos de `R/` (incluyendo los titulos de los 5
  datasets). En las cadenas de `R/` los caracteres acentuados se escriben
  con secuencias `\uXXXX`, de modo que el codigo R es ASCII puro (requisito
  de CRAN); la salida en pantalla conserva los acentos. Pendiente: correr
  `devtools::document()` para regenerar `man/*.Rd`.

## Diseño y presentacion

* Los 8 metodos `plot.*` ahora aceptan un argumento `pausar` (por defecto
  `interactive()`) y devuelven la lista de graficas ggplot en vez de
  `invisible(x)`, para poder componerlas con `patchwork`/`gridExtra`,
  guardarlas con `ggsave()` o usarlas en R Markdown sin bloquear con
  `readline()`.
* Las leyendas "INTERPRETACIÓN:" (texto fijo, identico sin importar el
  resultado) se renombraron a "CÓMO LEER ESTA GRÁFICA:" en las 74 graficas
  de diagnostico; ademas, la caption del VIF en `greenreg_rlm()` ahora reporta el
  VIF maximo real y la variable correspondiente en vez de un texto
  generico.
* Se eliminaron comentarios de codigo dirigidos al desarrollador ("Aqui
  unificamos el semaforo que creaste arriba", los 48 "CORRECCION N:", etc.)
  y se limpiaron 574 lineas con espacios en blanco al final, presentes en
  el paquete publicado en CRAN.
* Se mantiene "OKEY" / "ALERTA" como etiqueta de veredicto en las 9
  funciones (decision de diseño del paquete); se tradujeron al español los encabezados de
  tablas propias del paquete ("Estimate"/"Std. Error"/"z value" ->
  "Estimación"/"Error Estándar"/"valor z") y el encabezado mixto
  "(Ratio Test)" -> "(Prueba de Razón de Verosimilitud)"; los veredictos que
  afirmaban la hipótesis nula como hecho ("Los residuos son normales",
  "El modelo es ruido blanco") ahora dicen explícitamente "no se rechaza"
  en vez de darla por probada.

# GREENREG 0.1.1

Version publicada en CRAN, revisada como parte del comite de tesis.
