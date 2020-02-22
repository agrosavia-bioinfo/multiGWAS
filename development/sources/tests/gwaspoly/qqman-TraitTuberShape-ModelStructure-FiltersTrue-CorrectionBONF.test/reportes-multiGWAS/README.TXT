Hola Manfred,

Este script crea los diferentes reportes nombrados en la función "createReports"

- Para ejecutarlo primero se debe instalar las librerias de R llamadas al inicio del script, 
  las cuales son:
	library (stringr)
	library (dplyr)
	library (qqman)
	library (VennDiagram)

- Para ejecutarse necesita los archivos de entradas que están en el directorio "in"

- Las salidas las deja en el directorio "out", que se crea automaticamente

- La forma de ejecutarse desde una terminal linux (no he probado en windows) puede ser:
  $ ./gwas-summary-markdown.R
  ó
  $ Rscritp gwas-summary-markdown.R

- La función main simplemente llama a la funcion "createReports", la cual la he
  tratado de realizar de forma que llame desde allí la función de escritura
  externa, ya sea para las tablas con "write.table", o para los gráficos con
  "pdf" y "dev.off".
  Así creo que sería más fácil llevarlo a markdown ya que los scripts están
  listos para llamarse cuando se esté generando el markdown y creo que solo
  sería quitar las directivas de escritura externa o cambiarlas por las que 
  necesita para markdown.

- En la funcion "createReports" he tratado de crear un esquema básico de lo que
  llevaría el documento, título texto sencillo, y 5 secciones con subtitulos (H1)
  y texto sencillo. Así:
  :
  - Titulo Principal
  YYYYYYYYY ... texto
  Tabla de contenido

  - Sección 1. XXXXX
    YYYYYYYYYY ... texto
  - Sección 2. XXXXX
    YYYYYYYYYY ... texto
  - Sección 3. XXXXX
    YYYYYYYYYY ... texto
  - Sección 4. XXXXX
    YYYYYYYYYY ... texto
  - Sección 5. XXXXX
    YYYYYYYYYY ... texto

- Cualquier inquietud, me escribes o me llamas (tel 3176753141)

