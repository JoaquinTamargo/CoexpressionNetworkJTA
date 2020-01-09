TAREA 3. REDES DE COEXPRESIÓN GÉNICA.
===============================================

* Alumno/a: Joaquín Tamargo Azpilicueta
* E-mail: joatamazp@alum.us.es
* Biología Molecular de Sistemas, Grado en Bioquímica (Universidad de Sevilla)
------

## Tabla de Contenidos

1. [Introducción](#introduccion)
2. [Construcción de la red de co-expresión génica](#red)
3. [Análisis de la topología de la red](#topologia)
4. [Análisis de clustering](#clustering)
5. [Descripción de los clústeres identificados](#identificados)
6. [Conclusiones](#conclusion)
7. [Referencias](#referencias)

## Información complementaria
1. [Script completo](./script.r)
--------

# Red de coexpresión de genes regulados por factores de transcripción de primera oleada en respuesta a temperaturas bajas en *Arabidopsis thaliana*

## 1. Introducción 
<a class="anchor" id="introduccion"></a>

### 1.1. Objetivo del estudio

  Las plantas de regiones templadas, incluyendo *Arabidopsis thaliana* muestran un aumento en la tolerancia a la congelación en respuesta a temperaturas bajas (1). Aún no existe un consenso en el número e identidad de los genes de respuesta al frío en Arabidopsis (2), pero se cree que intervienen una serie de proteínas y moléculas crioprotectoras reguladas bajo genes conocidos como genes regulados por el frío (COR) (3). La vía más conocida de regulación de expresión de algunos de estos genes COR es la respuesta mediada por proteínas CBF (C-repeat binding factor), donde CBF1, CBF2 y CBF3 aumentan su expresión minutos después de exponer plantas de Arabidopsis a temperaturas bajas (3, 4). En este contexto, la expresión constitutiva de CBF1, CBF2 o CBF3 da lugar a una expresión alterada de unos 100 genes regulados por el frío que las hace tan resistentes aun sin exponerlas previamente a temperaturas bajas como las plantas silvestres expuestas a temperatura baja (1,6). 

Vogel et al. (2005) (3) describe que la expresión constitutiva de ZAT12, un factor de transcripción inducido con el frío paralelamente a CBF1, CBF2 y CBF3,  daba lugar a la inducción y represión de varios genes regulados por frío, de los cuales la mayoría también formaban parte del regulon de CBF2. De la misma forma que con CBF, la sobreexpresión de ZAT12 incrementó la tolerancia a temperaturas de congelación sin exponer previamente las plantas a temperaturas bajas, aunque el efecto no fuera tan intenso como con la sobreexpresión de CBF2. Esto sembró la duda de si el regulón de CBF también está regulado por otros factores de transcripción de "primera oleada" codificados por genes rápidamente inducidos en respuesta a temperaturas bajas paralelamente a CBF1, CBF2 y CBF3 y que todavía no hubieran sido descubiertos.

En el artículo, se identifican las proteínas inducidas por temperaturas bajas de no congelación y se analiza la regulación de estos por diferentes factores de transcripción de manera que pudo elaborar una red de coexpresión de genes inducidos por el frío. En el estudio de esta tarea, no obstante, el **objetivo** será, con los mismos resultados de transcriptómica obtenidos por los autores del artículo en revisión, realizar un análisis extensivo y exhaustivo del conjunto total de genes diferencialmente expresados con el fin de elaborar una red de coexpresión en la que participan diferentes factores de transcripción de primera oleada inducidos por el frío. Mediante técnicas de *clustering*, se obtienen de la red grupos de genes a los que se les hace un enriquecimiento de términos de ontología génica con el fin de conocer funciones potenciales de esos clústeres de genes.


### 1.2. Conceptos básicos en redes de co-expresión génica

El resultado de nuestro estudio será un grafo o red donde habrá una serie de nodos o vértices (que corresponden a las entidades de la red) unidos entre sí mediante aristas. Estas redes vienen definidas por parámetros locales y globales. En cuanto a los parámetros locales, el número de vecinos (el número de nodos conectados al nodo de estudio) se define también como el **grado de los nodos**. A partir de él se pueden calcular los valores de **transitividad**, centralidad, intermediación o excentricidad. De todos ellos, quizás el de mayor relevancia es el concepto de transitividad. Este (también llamado coeficiente de agrupamiento), representa el grado de agrupamiento en torno a un nodo y varía entre 0 (cuando todos los vecinos no son vecinos entre sí) y 1 (cuando todos los vecinos del nodo de interés son vecinos entre sí). Viene definido por el grado de los nodos (d) de un nodo v y por las aristas entre los vecinos de v (N):

<font size="1">.</font><CENTER>$t_v=\frac{N_v}{\frac{d_v(d_v-1)}{2}}$</CENTER>

Por otro lado, uno de los parámetros globales más importantes es la **distribución del grado de los nodos**. Esta se define como la probabilidad de encontrar un nodo de grado k de forma aleatoria en la red. Esta propiedad de la red es importante en tanto que la distribución del grado de los nodos en las *redes libres de escala* (que son la mayoría de redes biológicas caracterizadas experimentalmente) siguen una distribución potencial negativa. Esto es:

<font size="1">.</font><CENTER>$P(k)=k*K^\alpha$   |  $ \alpha<0$</CENTER>

Esto implica que encontrar nodos con muchas conectividades es poco probable, mientras que es mucho más probable encontrar nodos con pocos vecinos. De esta forma, en estas redes hay nodos concentradores de muuchas aristas, conocidos como *hubs*. Esto permite que en estas redes tenga lugar una propagación rápida de la información al tener diámetros, radios y longitudes medias de caminos mínimos pequeñas si se comparan con redes aleatorias. Así, las redes libres de escala son robustas a ataques aleatorios, pero frágiles a ataques dirigidos (hacia los hubs).

Una *red de mundo pequeño* se define como una red libre de escala que presenta un alto coeficiente de agrupamiento, y donde el camino entre dos nodos cualesquiera es pequeño.

Se puede ampliar la información aquí sintetizada visitando [esta página](https://frannetworks.shinyapps.io/red_trenes/) (7).

### 1.3. Diseño experimental

Para la **identificación de los factores de transcripción inducidos por el frío**, en el artículo se usaron los datos de los estudios de expresión por microarrays de Kilian *et al.* (2007) (8), disponibles en [GSE5620](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse5620) y [GSE5621](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse5621)). Brevemente, ese estudio consistió en un experimento de 24h de duración donde se tomaron los brotes de plantas de *Arabidopsis thaliana Wild Type (col-0)* expuestas (o no) al frío (4 ºC). Con esto, los autores del artículo en discusión fueron capaces de identificar 2397 genes expresados diferencialmente (con un fold-change, FC = 4 y una tasa de descubrimientos falsos, FDR = 0,01). De estos, se seleccionaron 1336, que correspondían a los genes inducidos por frío: 174 codificaban para factores de transcripción, de los que 27 genes se inducían en paralelo a CBF1, CBF2 y CBF3. 17 de estos 27 genes se indujeron rápidamente bajo las condiciones experimentales. Se consiguieron líneas transgénicas que sobreexpresaban bajo el promotor del RNA 35S del virus del mosaico de la coliflor 11 de estos factores de transcripción: MYB73, CRF2, RAV1, CRF3, ERF5, DEAR1, MYB44, CZF1, ZAT10, ZF y HSFC1 (Figura 1). No se describieron ensayos para probar los efectos de la inserción del promotor CaMV 35S. Este promotor CaMV 35S es particularmente controvertido puesto que algunos fragmentos largos de este promotor contienen un marco de lectura abierto (ORF) codificante para la proteína multifuncional P6 que puede desencadenar efectos fenotípicos no deseados (9).

<img src="./00_Images/Imagen1.png" alt="figura1" height="720" width="360">

<span style="font-size:0.85em">**Figura 1.** Diagrama de la selección de los genes para el estudio.</span>

Para este estudio, se hicieron dos réplicas independientes para cada una de las condiciones. Se usaron como **control**  plantas sometidas a las mismas condiciones que las **plantas modificadas** que sobreexpresaban los 11 genes anteriormente mencionados con el fin de determinar la interrelación entre la red de transcritos cuando se expresan constitutivamente los 11 distintos factores de transcripción. Las semillas de las plantas silvestres (WT) y de la generación T3 de plantas transgénicas fueron germinadas a 22 ºC bajo condiciones de esterilidad en medio sólido de agar y Phytoblend 0,8% con medio de cultivo de Gamborg B5 sin sacarosa durante 12 días a luz constante de ~100 µmol m$^–$$^2$ s$^–$$^1$. Se extrajo el RNA total de las plántulas y se biotinizó usando el kit de etiquetado de Affymetrix para IVT. Se hibridó con el Affymetrix Arabidopsis ATH1 GenChip (los datos de las anotaciones se pueden encontrar en la página de Bioconductor de [ath1121501.db](https://bioconductor.org/packages/release/data/annotation/html/ath1121501.db.html)) (20). Los datos crudos fueron depositados en en *Gene Expression Omnibus* ([GEO](https://www.ncbi.nlm.nih.gov/geo/)) del NCBI bajo el número de identificación [GSE55907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55907). 

En los datos crudos también se incluyen los datos de los análisis de microarrays de una línea WT de control (incubado a 22 ºC) y de una línea WT incubada a 4 ºC, con el fin de identificar genes diferencialmente expresados con temperaturas bajas y seleccionar estos genes dentro de los activados o reprimidos cuando se expresa constitutivamente cualquiera de los factores de transcripción de las líneas transgénicas. En el estudio que se presenta en esta tarea, se omitirá esta selección porque se hará el estudio a nivel global y no a nivel únicamente de genes activados en respuesta al frío, y por tanto no se tendrán en consideración los datos correspondientes a GSM1348266, GSM1348267, GSM1348268 y GSM1348269.

Este estudio se busca la relación potencial en los genes modulados por 11 factores de transcripción de primera oleada de respuesta a frío. No obstante, la red de regulación de aclimatación al frío es muy compleja, y solamente 11 factores de transcripción de los 27 factores de transcripción de la primera oleada identificados han sido probados en este artículo, por lo que es muy posible que existan factores de transcripción diferentes a los CBF que participen en la regulación de la expresión de genes protectores del frío (10) y no se han tenido en consideración.

### 1.4. Flujo de trabajo

Para la revisión, se descargarán y se analizarán los datos del estudio usando [RStudio](https://rstudio.com/) (21). Se leerán mediante el paquete [*affy*](https://www.bioconductor.org/packages/release/bioc/html/affy.html) (22) los datos crudos del análisis de microarrays correspondientes a la identificación [GSE55907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55907). Se analizará la calidad de los datos mediante el control de calidad implementado en [*simpleaffy*](https://bioconductor.org/packages/release/bioc/html/simpleaffy.html) (23), se normalizarán y se extraerán los niveles de expresión para elaborar una lista con los genes diferencialmente expresados mediante un contraste donde se comparan los niveles de expresión en el control (WT) y los que sobreexpresan algún factor de transcripción, usando el paquete [*limma*](https://www.bioconductor.org/packages/release/bioc/html/limma.html) (24). Se determina la coexpresión a partir de un análisis de la correlación entre los perfiles de expresión génica. A continuación, se elige un umbral de corte en la correlación para generar la matriz de adyacencia con la que se elaborará la red. Para ello, se generan las redes para diferentes umbrales de corte y se estudian los valores de conectividad media de los nodos así como su ajuste a una potencial negativa (lo que es indicadora de que las redes son libres de escala), usando para ello el paquete [*igraph*](https://cran.r-project.org/web/packages/igraph/index.html) (25). Habiendo escogido el umbral y generado la red, se visualizará mediante el programa de visualización de redes [Cytoscape](https://cytoscape.org/) (26). Se analizará la topología de la red para: 1) determinar la propiedad de red libre de escala, mediante el test de Kolmogorov-Smirnov y mediante una regresión lineal y 2) la determinación de red de mundo pequeño. Se identificarán y caracterizarán los *hubs* y se hará un clustering de la red. Para ello, se comparará la bondad de los métodos jerárquico y de partición en torno a centroides y se escogerá el de mejor resultado. De los clústeres obtenidos, así como de los hubs de la red, se hará un enriquecimiento de términos de ontología de genes mediante [GOrilla](http://cbl-gorilla.cs.technion.ac.il/) (27, 28) y [ReViGO](http://revigo.irb.hr/) (29).

<img src="./00_Images/jobflow.png" alt="figura1" height="720" width="1080">

## 2. Construcción de la red de co-expresión génica 
<a class="anchor" id="red"></a>

### 2.1. Procesamiento de los datos y estimación de los niveles de expresión
En primer lugar, se descargaron los datos de [GSE55907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55907) y se eliminaron los ficheros correspondientes a las dos líneas WT usadas para los análisis de expresión génica diferencial en condiciones de frío (GSM1348266, GSM1348267, GSM1348268 y GSM1348269). Se leyeron con *ReadAffy* (22) y se efectuó un control de calidad (véase el script) con el algoritmo implementado en el paquete *simpleaffy* (23) antes de normalizar los valores de fluorescencia de fondo y extraer los niveles de expresión (mediante *rma*, que los extrae en log$_2$). En el control de calidad se vio una advertencia para los valores del porcentaje del total de sondas que presentaban fluorescencia, pero estos valores eran muy parecidos en todas las muestras y ninguna muestra presentaba advertencias para las sondas de control de degradación. Se asumió, por tanto, un buen estado general de las muestras.


```R
suppressMessages(library(affy))

microarray.raw.data <- ReadAffy(verbose = TRUE)

microarray.preprocessed.data <- rma(microarray.raw.data)
expression.levels <- exprs(microarray.preprocessed.data)
```

    1 reading /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348270_Col0_1.CEL.gz ...instantiating an AffyBatch (intensity a 506944x24 matrix)...done.
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348270_Col0_1.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348271_Col0_2.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348272_CRF2_1.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348273_CRF2_2.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348274_CRF3_1.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348275_CRF3_2.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348276_CZF1_1.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348277_CZF1_2.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348278_DEAR1_1.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348279_DEAR1_2.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348280_ERF5_1.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348281_ERF5_2.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348282_HSFC1_1.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348283_HSFC1_2.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348284_MYB44_1.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348285_MYB44_2.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348286_MYB73_1.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348287_MYB73_2.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348288_RAV1_1.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348289_RAV1_2.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348290_ZAT10_1.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348291_ZAT10_2.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348292_ZF_1.CEL.gz
    Reading in : /Users/ParaisoNoFiscal/Biologia Molecular de Sistemas - Jupyter/T3-Tamargo_Azpilicueta_Joaquin/GSM1348293_ZF_2.CEL.gz


    Warning message:
    “replacing previous import ‘AnnotationDbi::tail’ by ‘utils::tail’ when loading ‘ath1121501cdf’”
    Warning message:
    “replacing previous import ‘AnnotationDbi::head’ by ‘utils::head’ when loading ‘ath1121501cdf’”
    
    


    Background correcting
    Normalizing
    Calculating Expression


En primer lugar, se previsualizó el efecto de las mutaciones con la elaboración de gráficos de dispersión. Para ello, y para la obtención de la matriz de correlación a partir de la cual se obtienen las distancias para la generación de la red, es necesario obtener una matriz de expresión media (mean.expression) que disponga en cada columna la expresión media de cada réplica para cada una de las condiciones.

Se harán comparaciones del control a 22 ºC con las 11 líneas transgénicas que expresan constitutivamente los factores de transcripción. Así, se podrá determinar el efecto general en el transcriptoma cuando se expresan los mencionados factores de transcripción de primera oleada de respuesta al frío.


```R
col0 <- (expression.levels[,1]+expression.levels[,2])/2
crf2 <- (expression.levels[,3]+expression.levels[,4])/2
crf3 <- (expression.levels[,5]+expression.levels[,6])/2
czf1 <- (expression.levels[,7]+expression.levels[,8])/2
dear <- (expression.levels[,9]+expression.levels[,10])/2
erf5 <- (expression.levels[,11]+expression.levels[,12])/2
hsfc1 <- (expression.levels[,13]+expression.levels[,14])/2
myb44 <- (expression.levels[,15]+expression.levels[,16])/2
myb73 <- (expression.levels[,17]+expression.levels[,18])/2
rav1 <- (expression.levels[,19]+expression.levels[,20])/2
zat10 <- (expression.levels[,21]+expression.levels[,22])/2
zf <- (expression.levels[,23]+expression.levels[,24])/2

mean.expression <- matrix(c(col0, crf2, crf3, czf1, dear, erf5, hsfc1, myb44, myb73, rav1, zat10, zf), ncol=12)
experiment.name <- c("col0",  "crf2", "crf3", "czf1", "dear", "erf5", "hsfc1", "myb44", "myb73", "rav1", "zat10", "zf")
colnames(mean.expression) <- experiment.name

par(mfrow=c(3,4))
plot(mean.expression[,"col0"], mean.expression[,"crf2"], pch = 20, col = "black", xlab = "WT", ylab = "crf2")
plot(mean.expression[,"col0"], mean.expression[,"crf3"], pch = 20, col = "black", xlab = "WT", ylab = "crf3")
plot(mean.expression[,"col0"], mean.expression[,"czf1"], pch = 20, col = "black", xlab = "WT", ylab = "czf1")
plot(mean.expression[,"col0"], mean.expression[,"dear"], pch = 20, col = "black", xlab = "WT", ylab = "dear")
plot(mean.expression[,"col0"], mean.expression[,"erf5"], pch = 20, col = "black", xlab = "WT", ylab = "erf5")
plot(mean.expression[,"col0"], mean.expression[,"hsfc1"], pch = 20, col = "black", xlab = "WT", ylab = "hsfc1")
plot(mean.expression[,"col0"], mean.expression[,"myb44"], pch = 20, col = "black", xlab = "WT", ylab = "myb44")
plot(mean.expression[,"col0"], mean.expression[,"myb73"], pch = 20, col = "black", xlab = "WT", ylab = "myb73")
plot(mean.expression[,"col0"], mean.expression[,"rav1"], pch = 20, col = "black", xlab = "WT", ylab = "rav1")
plot(mean.expression[,"col0"], mean.expression[,"zat10"], pch = 20, col = "black", xlab = "WT", ylab = "zat10")
plot(mean.expression[,"col0"], mean.expression[,"zf"], pch = 20, col = "black", xlab = "WT", ylab = "zf")
```


![png](output_6_0.png)


El número de genes diferencialmente expresados difirió enormemente entre los 11 factores de transcripción constitutivamente expresados. Estos resultados son coherentes con el artículo en revisión, donde se describe un gran efecto de HSFC1 en la transcripción de diferentes genes. Los autores señalan que CRF2, CRF3, DEAR1, ERF5, MYB44, MYB73 y RAV1 no afectaron a la expresión de ninguno de los genes regulados por el frío (COR). En esta previsualización se aprecia un efecto mucho menor sobre la expresión génica en estos genotipos que en el de HSFC1 o CZF1, por ejemplo. Aun así se mantendrán los datos con la intención de hacer un estudio de los efectos a nivel global de estos factores de transcripción, y para conocer si participan en procesos similares las proteínas que potencialmente coexpresan. Para poder determinar este efecto, y conocer en profundidad qué genes están involucrados en la respuesta, se debe hacer una selección de los genes expresados de forma diferencial (DEGs).

### 2.2. Análisis de expresión génica diferencial

[LIMMA](https://bioconductor.org/packages/release/bioc/html/limma.html) (*LInear Models for Microarray Analysis*) (24) es un paquete de Bioconductor para el análisis de análisis de datos de expresión génica de microarrays, concretamente para la evaluación de la expresión diferencial. Para poder proceder a la selección de genes diferencialmente expresados, primero se genera una matriz que represente el diseño experimental, asignando a cada condición experimental un número entero. El control recibe el número 1 y las líneas transgénicas reciben el 2 (CRF2), 3 (CRF3), 4 (CZF1), 5 (DEAR), 6 (ERF5), 7 (HSFC), 8 (MYB4), 9 (MYB7), 10 (RAV1), 11 (ZAT10) y 12 (ZF).

Se asignará uno de estos 12 números a cada una de las muestras con su etiqueta correspondiente, usando la función *model.matrix* contenida en el paquete LIMMA . A continuación, atendiendo al diseño experimental, con la función *lmFit*, se ajusta la estimación de los niveles de expresión de cada gen a un modelo lineal. Básicamente, esto permitirá obtener la media de las réplicas para cada condición. 

Para determinar los genes que se expresan diferencialmente, se comparan los niveles de expresión del único control con los de las líneas transgénicas. Para especificar los contrastes se usa la función *makeContrasts*, que construye la matriz de contrastes del conjunto de datos que se especifican. La función *eBayes* permite ordenar los genes en función de su expresión diferencial.

Existen diferentes criterios para definir la **expresión génica diferencial**. En esta tarea, se usará el método mixto de fold-change y p-valor. El método del fold-change es adecuado al tratarse de organismos modelo y disponer de un limitado número de réplicas para cada experimento independiente, pero este método no tiene controles de falsos positivos, por lo que adicionalmente se tomará un umbral de significancia estadística para la expresión. Los genes diferencialmente expresados serán estimados con la función DEG.selection.

La función DEG.selection recibe como entrada los datos procesados (processed.data), el número de comparaciones (number.comparisons), el número de genes (number.genes), el umbral de FC (fold.change.threshold) y el del p-valor (log.pvalue.threshold), ambos en log$_{10}$. Con esta función, se estimarán los genes diferencialmente expresados a partir de contrast.results para las 11 comparaciones, con 22810 genes de la placa de microarrays, un umbral de FC de 2 (log$_2$(FC) = 1) y un p-valor de 0,05 (log$_{10}$(p-value) = 1.3). 

Las funciones que se usan a continuación para calcular la correlación en la expresión de los DEGS requieren que en las filas se encuentren las condiciones y en las columnas aparezcan los genes de cada sonda, al contrario que como aparecían en la matriz de expresión media. Es entonces necesario obtener la matriz traspuesta a mean.expression y seleccionar en ella sólo los genes correspondientes a los genes diferencialmente expresados.


```R
suppressMessages(library(limma))

experimental.design <- model.matrix(~ -1+factor(sort(c(seq(from=1, to=12),seq(from=1, to=12)))))
colnames(experimental.design) <- experiment.name

linear.fit <- lmFit(expression.levels, experimental.design)

contrast.matrix <- makeContrasts(crf2-col0, crf3-col0, czf1-col0, dear-col0, erf5-col0,
                                 hsfc1-col0, myb44-col0, myb73-col0, rav1-col0,
                                 zat10-col0, zf-col0, 
                                 levels= c("col0", "crf2","crf3","czf1","dear",
                                           "erf5","hsfc1","myb44","myb73","rav1",
                                           "zat10","zf"))

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)


DEG.selection <- function(processed.data,number.comparisons,number.genes,
                          fold.change.threshold,log.pvalue.threshold)
{
  
  ## Inicializar un vector vacío de caracteres
  DEGs <- character()
  
  ## Bucle for para recorrer todas las comparaciones
  for(i in 1:number.comparisons)
  {
    ## Generación de una tabla con información de expresión diferencial 
    DEGs.partial <- topTable(processed.data, number=number.genes,coef=i)
    
    ## Extracción del FC, del exponente del q-valor y del nombre de las sondas
    ## de los genes diferencialmente expresados
    
    fold.change <- DEGs.partial[["logFC"]]
    log.p.value <- -log10(DEGs.partial[["adj.P.Val"]])
    probe.names <- rownames(DEGs.partial)
    
    ## Los genes activados se determinan como aquellos que superan el umbral del
    ## FC con la significancia estadística especificada
    
    activated <- (fold.change > fold.change.threshold) & (log.p.value > log.pvalue.threshold)

    ## Los genes reprimidos se determinan como aquellos que quedan por debajo del 
    ## umbral de FC con la significancia estadística especificada 
    
    inhibited <- (fold.change < - fold.change.threshold) & (log.p.value > log.pvalue.threshold)
    
    ## Se extraen los nombres de las sondas para los genes activados e inhibidos
    
    activated.genes <- probe.names[activated]
    inhibited.genes <- probe.names[inhibited]
    
    ## Se añaden los genes activados e inhibidos en la comparación actual al acumulador de
    ## DEGs
    
    DEGs <- c(DEGs, activated.genes, inhibited.genes)
  }
  
  ## La función tiene como salida el vector acumulador de DEGs después de eliminar las
  ## repeticiones
  
  return(unique(DEGs))
}

DEG.selection.result <- DEG.selection(processed.data = contrast.results,
                               number.comparisons = 11, number.genes = 22810, 
                               fold.change.threshold = 1, log.pvalue.threshold = 1.3)

traspuesta<-t(mean.expression)
rownames(traspuesta) <- experiment.name
colnames(traspuesta) <- rownames(expression.levels)
diff.expr <- traspuesta[,DEG.selection.result]
dim(diff.expr)
```


<ol class=list-inline>
	<li>12</li>
	<li>2488</li>
</ol>



### 2.3. Análisis de correlación entre los perfiles de expresión de los genes expresados de forma diferencial

El criterio que se seguirá para determinar la coexpresión de dos genes se basará en la correlación entre los perfiles de expresión (correspondientes a los distintos experimentos) para los genes expresados de forma diferencial. De esta forma, se calculará la matriz de correlaciones génicas a partir de la matriz diff.expr usando la función *cor*. Con esta función, se obtiene una matriz cuadrada con tantas columnas y filas como DEGs. Será con esta matriz con la que elaboremos la matriz de distancias para los métodos de clustering.


```R
correlation.matrix <- cor(diff.expr)
dim(correlation.matrix)
```


<ol class=list-inline>
	<li>2488</li>
	<li>2488</li>
</ol>



### 2.4. Construcción de la red de coexpresión: selección de la matriz de adyacencia

Es necesario escoger un umbral de corte en la correlación para poder determinar si se encuentran co-expresados de forma significativa. Para ello, se genera un vector que almacena unos umbrales que van a ser recorridos con un bucle for. Para cada valor del umbral, se obtendrá la matriz de adyacencia y se elaborará con ella una red a la que se le calculará la conectividad media de los nodos. 

La mayoría de redes biológicas caracterizadas de forma experimental cumplen con la propiedad de ser redes libre de escala; es decir, que su distribución del grado de los nodos sigue una ley potencial negativa. Por ello, se calcularán para cada umbral los valores de R2 para el ajuste lineal sobre la transformada logarítmica de la distribución del grado de los nodos. Un R2 mayor indicará un mejor ajuste a la distribución lineal de la transformada logarítmica, lo que indicará que sigue una ley potencial, y será más probable que se trate de una red libre de escala. Este punto se desarrollará con algo más de detalle en el punto 3.1.2.

En conclusión, el umbral debe escogerse de manera que se obtengan valores de R2 altos y valores de conectividad media (al menos) mayores que 10.


```R
suppressMessages(library(igraph))

## Crear un vector con los umbrales y vectores vacíos para almacenar las
## conectividades medias y los R2.

thresholds <- seq(from=0.70,to=0.99,by=0.01)
mean.connectivities <- vector(length=length(thresholds))
scale.free.R2 <- vector(length=length(thresholds))
gene.correlation <- correlation.matrix

par(mfrow=c(1,1))
for(i in 1:length(thresholds))
{
  # print(thresholds[i])
  
  ## Construcción de la red correspondiente al umbral actual
  
  ## MATRIZ DE ADYACENCIA
  current.adjacency <- (gene.correlation > thresholds[i] & gene.correlation < 1)
  
  ## RED
  threshold.network <- graph.adjacency(current.adjacency, mode="undirected")
  
  ## GRADO DE LOS NODOS
  node.degrees <- degree(threshold.network)
  
  ## CONECTIVIDAD MEDIA (Calculada como la media del grado de los nodos)
  mean.connectivities[i] <- mean(node.degrees)
  
  ##
  ## COMPROBACIÓN DE QUE SE TRATA DE UNA RED LIBRE DE ESCALA
  ## Comprobación gráfica (se grafica el grado de los nodos)
  # h <- hist(node.degrees)
  
  ## Frecuencia absoluta del grado de los nodos
  degree.frequencies <- table(node.degrees)
  
  ## Si sigue una exponencial negativa, entonces y=b*k^a
  ## De esta forma, tomando logaritmos: log(y)=a*log(k)+log(b)
  
  ## Por tanto, se hace una determinación de una regresión lineal
  ## para las frecuencias de los grados en escala logarítmicas, 
  ## extrayendo el primer valor (que es 0).
  ##

  lm.r <- lm(log10(as.numeric(names(degree.frequencies[-1]))) ~ log10(degree.frequencies[-1]))
  
  ## Se extrae el R cuadrado como una medida de ajuste a la propiedad
  ## libre de escala
  s.lm <- summary(lm.r)
  scale.free.R2[i] <- s.lm[["adj.r.squared"]]
}

## Se grafican las conectividades medias y el R2 calculado en los
## pasos anteriores en función de los distintos umbrales almacenados
## en la variable "thresholds". A continuación, se muestran los valores.

par(mfrow = c(1,2))
plot(thresholds,mean.connectivities,type="o",col="red",lwd=3,xlab="Correlation Threshold",ylab="Mean connectivity")
plot(ylim=c(0,1),thresholds,scale.free.R2,type="o",col="blue",lwd=3,xlim=c(0.70,0.99),xlab="Correlation Threshold",ylab="Scale Free Model Fit (R²)")

names(mean.connectivities) <- thresholds
names(scale.free.R2) <- thresholds

mean.connectivities
scale.free.R2
```


<dl class=dl-horizontal>
	<dt>0.7</dt>
		<dd>215.309485530547</dd>
	<dt>0.71</dt>
		<dd>203.228295819936</dd>
	<dt>0.72</dt>
		<dd>191.327974276527</dd>
	<dt>0.73</dt>
		<dd>179.838424437299</dd>
	<dt>0.74</dt>
		<dd>168.651125401929</dd>
	<dt>0.75</dt>
		<dd>157.638263665595</dd>
	<dt>0.76</dt>
		<dd>147.006430868167</dd>
	<dt>0.77</dt>
		<dd>136.577974276527</dd>
	<dt>0.78</dt>
		<dd>126.20498392283</dd>
	<dt>0.79</dt>
		<dd>116.324758842444</dd>
	<dt>0.8</dt>
		<dd>106.725884244373</dd>
	<dt>0.81</dt>
		<dd>97.4646302250804</dd>
	<dt>0.82</dt>
		<dd>88.435691318328</dd>
	<dt>0.83</dt>
		<dd>79.8392282958199</dd>
	<dt>0.84</dt>
		<dd>71.467845659164</dd>
	<dt>0.85</dt>
		<dd>63.5233118971061</dd>
	<dt>0.86</dt>
		<dd>55.9469453376206</dd>
	<dt>0.87</dt>
		<dd>48.6864951768489</dd>
	<dt>0.88</dt>
		<dd>41.717845659164</dd>
	<dt>0.89</dt>
		<dd>35.2942122186495</dd>
	<dt>0.9</dt>
		<dd>29.1109324758842</dd>
	<dt>0.91</dt>
		<dd>23.4453376205788</dd>
	<dt>0.92</dt>
		<dd>18.3336012861736</dd>
	<dt>0.93</dt>
		<dd>13.7596463022508</dd>
	<dt>0.94</dt>
		<dd>9.85691318327974</dd>
	<dt>0.95</dt>
		<dd>6.55787781350482</dd>
	<dt>0.96</dt>
		<dd>3.87540192926045</dd>
	<dt>0.97</dt>
		<dd>1.89951768488746</dd>
	<dt>0.98</dt>
		<dd>0.688906752411576</dd>
	<dt>0.99</dt>
		<dd>0.104501607717042</dd>
</dl>




<dl class=dl-horizontal>
	<dt>0.7</dt>
		<dd>0.0723085066741898</dd>
	<dt>0.71</dt>
		<dd>0.0991799463175921</dd>
	<dt>0.72</dt>
		<dd>0.123550209753733</dd>
	<dt>0.73</dt>
		<dd>0.137686004733942</dd>
	<dt>0.74</dt>
		<dd>0.187300481641568</dd>
	<dt>0.75</dt>
		<dd>0.227895524315779</dd>
	<dt>0.76</dt>
		<dd>0.278735269851777</dd>
	<dt>0.77</dt>
		<dd>0.307988246970737</dd>
	<dt>0.78</dt>
		<dd>0.370512347681642</dd>
	<dt>0.79</dt>
		<dd>0.426008593533115</dd>
	<dt>0.8</dt>
		<dd>0.496028308841082</dd>
	<dt>0.81</dt>
		<dd>0.550026887386798</dd>
	<dt>0.82</dt>
		<dd>0.613962705384679</dd>
	<dt>0.83</dt>
		<dd>0.658894297130283</dd>
	<dt>0.84</dt>
		<dd>0.703781439547664</dd>
	<dt>0.85</dt>
		<dd>0.704751900620886</dd>
	<dt>0.86</dt>
		<dd>0.716353063976917</dd>
	<dt>0.87</dt>
		<dd>0.727562951915442</dd>
	<dt>0.88</dt>
		<dd>0.74030114993558</dd>
	<dt>0.89</dt>
		<dd>0.726044113656103</dd>
	<dt>0.9</dt>
		<dd>0.730380624148255</dd>
	<dt>0.91</dt>
		<dd>0.675080108227159</dd>
	<dt>0.92</dt>
		<dd>0.727310203275407</dd>
	<dt>0.93</dt>
		<dd>0.705275686964636</dd>
	<dt>0.94</dt>
		<dd>0.684048425469639</dd>
	<dt>0.95</dt>
		<dd>0.675698434273759</dd>
	<dt>0.96</dt>
		<dd>0.597840284375962</dd>
	<dt>0.97</dt>
		<dd>0.674416116161225</dd>
	<dt>0.98</dt>
		<dd>0.85594878520393</dd>
	<dt>0.99</dt>
		<dd>0.767770698272887</dd>
</dl>




![png](output_12_2.png)


Para un valor del umbral 0.88, se tiene un R$^2$ relativamente alto si se compara con el resto de valores (R$^2$=0.74) y una conectividad de más de 41 vecinos de media. De esta manera, 0.88 se toma como el umbral para la elaboración de la matriz de adyacencia, a partir de la cual se genera la red de coexpresión génica, que se puede visualizar en Cytoscape.


```R
adjacency.88 <- (gene.correlation > 0.88) & (gene.correlation < 1)
gene.coexpression.network.88 <- graph.adjacency(adjacency.88, mode="undirected")

write.graph(gene.coexpression.network.88,file="ath_gene_coexpression_network_88.gml",format="gml")
```

Se muestra a continuación el grafo que resultó del procesamiento anterior. Se señalan sobre él algunas proteínas de las que se habla a lo largo del texto.

<img src="./00_Images/jta_selected.png" alt="Smiley face" height="720" width="1080">

## 3. Análisis de la topología de la red
<a class="anchor" id="topologia"></a>

### 3.1. Determinación de propiedad de red libre de escala
#### 3.1.1. Método de Kolmogorov-Smirnov

Para conocer si una red es libre de escala se estudia si la distribución del grado de los nodos sigue una exponencial negativa (véase apartado 1.2). Para ello, con el test de Kolmogorov-Smirnov se hace un contraste de hipótesis donde: la hipótesis nula (H0): la distribución del grado de los nodos sigue una potencial y la hipótesis alternativa (H1): la distribución de los grados de los nodos no sigue una distribución potencial.

Para este contraste, se extrae la distribución del grado de los nodos mediante la función *degree.distribution* y se hace la comprobación de la presencia de la ley de potencia mediante *power.law.fit*, extrayéndose el valor correspondiente al método de Kolmogorov-Smirnov (KS.p):


```R
network.degree.distribution <- degree.distribution(gene.coexpression.network.88)
fit.scale.free <- power.law.fit(network.degree.distribution)
fit.scale.free[["KS.p"]]
```


0.948718357566543


Como 0.949 es mucho mayor que el nivel de significación, no se puede descartar la hipótesis nula: la distribución
del grado de los nodos sigue una potencial negativa y se trata por tanto de una red libre de escala.

#### 3.1.2. Método de regresión lineal

Las transformaciones logarítmicas transforman las potenciales negativas en rectas con pendiente negativa:

<font size="1">.</font><CENTER>$P(k)=\beta*K^\alpha\implies log_{10}(P(k))=\alpha*log_{10}(K)+log_{10}(\beta)$</CENTER>

Por tanto, para comprobar el ajuste a una ley potencial negativa sería suciente realizar una regresión lineal sobre los datos transformados y comprobar que la R$^2$ (el porcentaje de la varianza explicada por el ajuste) es alta y el p-valor es suficientemente bajo. De esta manera, se obtiene la frecuencia absoluta del grado de los nodos. Para poder aplicar logaritmos, es necesario eliminar los ceros. Se hace una regresión lineal con la función *lm*.




```R
par(mfrow=c(1,1))
network.degrees <- degree(gene.coexpression.network.88)

degree.frequencies <- table(network.degrees)
degree.frequencies.no.0 <- degree.frequencies[-1]

log10.degrees.frequencies <- log10(degree.frequencies.no.0)
log10.node.degrees <- log10(as.numeric(names(degree.frequencies.no.0)))

lm.r <- lm(log10.degrees.frequencies ~ log10.node.degrees)
summary(lm.r)

degree.histogram <- hist(network.degrees,
                         freq=FALSE, col="light blue",
                         xlab="Node Degrees",
                         ylab="Probability", main="Distribución del grado de los nodos")
```


    
    Call:
    lm(formula = log10.degrees.frequencies ~ log10.node.degrees)
    
    Residuals:
         Min       1Q   Median       3Q      Max 
    -0.73527 -0.17691  0.04115  0.18332  0.65997 
    
    Coefficients:
                       Estimate Std. Error t value Pr(>|t|)    
    (Intercept)         2.58821    0.07816   33.11   <2e-16 ***
    log10.node.degrees -1.01115    0.03871  -26.12   <2e-16 ***
    ---
    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    
    Residual standard error: 0.2578 on 238 degrees of freedom
    Multiple R-squared:  0.7414,	Adjusted R-squared:  0.7403 
    F-statistic: 682.3 on 1 and 238 DF,  p-value: < 2.2e-16




![png](output_18_1.png)


Se obtiene que para el valor ajustado de R cuadrado que escogimos como umbral (0.7403), el p-valor de la regresión lineal de la transformación logarítmica es mucho menor que 0.001, y que el valor para $\alpha$ es -1.01. Podemos ver el ajuste a una recta de pendiente negativa en el siguiente gráfico. En conclusión, la red escogida es una red libre de escala. Faltará conocer si esta red libre de escala es de mundo pequeño.


```R
lm.r <- lm(log10.degrees.frequencies ~ log10.node.degrees)
plot(x = log10.node.degrees, y = log10.degrees.frequencies, pch=20, col = "grey", ylim=c(0,3),
     main="Regresión de la frecuencia de los grados de los nodos 
     con transformación logarítmica", xlab="log(Grado de los nodos)", ylab="log(Frecuencia de grados)")
abline(lm.r, col="blue")
```


![png](output_20_0.png)


### 3.2. Determinación de propiedad de red de mundo pequeño

La condición que debe cumplir una red de mundo pequeño es que sus nodos, aunque sean muy abundantes, estén muy conectados y sea posible hacer caminos cortos para llegar de un nodo a otro cualquiera de la red. Los valores de la transitividad (o coeficiente de agrupamiento) y la longitud media del camino mínimo son principalmente los indicadores que permitirán hacer la determinación de que la red sea de mundo pequeño.

La transitividad representa el grado de agrupamiento en torno a un nodo, y toma valores entre 0 (todos los vecinos del nodo de interés no son vecinos entre sí) y 1 (todos los nodos vecinos del nodo de interés son vecinos entre sí) y puede ser estimada con la función *transitivity* a nivel global.

La longitud media del camino mínimo puede estimarse con *average.path.length*, y corresponde a la media de la distancia mínima que debe recorrerse para alcanzar de un nodo a otro cualquiera de la red. Se entiende que cuanto menor es este valor, es más probable que la red sea de mundo pequeño.


```R
network.clustering.coefficient <- transitivity(gene.coexpression.network.88,type="global")
network.clustering.coefficient
average.path.length(gene.coexpression.network.88, directed=FALSE)
```


0.718225641836312



6.46080827131269


Para comprobar que la generación de la red no ha tenido lugar de forma puramente aleatoria, se estimará la probabilidad de obtener una red similar pero con longitudes medias del camino mínimo inferior a la estudiada mediante el método de Albert-Lászlo  Barabasi. Según este, se parte de una pequeña red donde todos los nodos están unidos entre sí y se van conectando nuevos nodos a los  ya existentes, pero no de forma arbitraria, sino con una preferencia mayor por los nodos de alto grado. Para este propósito se usará la función *barabasi.game*, que recibe como entrada el tamaño de la red a crear (el número de nodos) y un vector numérico que da el número de aristas a añadir en cada paso, y devuelve una red libre de escala.

En primer lugar, se evalúa la red de coexpresión génica generada para comprobar el número de nodos y aristas de la red.


```R
gene.coexpression.network.88
```


    IGRAPH b50a436 UN-- 2488 51897 -- 
    + attr: name (v/c)
    + edges from b50a436 (vertex names):
     [1] 248253_at--260555_at   248253_at--260132_s_at 248253_at--253996_at  
     [4] 248253_at--247739_at   248253_at--258623_at   248253_at--251707_at  
     [7] 248253_at--255749_at   248253_at--265899_s_at 248253_at--257058_at  
    [10] 248253_at--251537_at   248253_at--246457_at   248253_at--259941_s_at
    [13] 248253_at--250495_at   248253_at--249253_at   248253_at--256016_at  
    [16] 260555_at--260132_s_at 260555_at--250179_at   260555_at--258821_at  
    [19] 260555_at--249975_s_at 260555_at--253202_at   260555_at--254429_at  
    [22] 260555_at--258623_at   260555_at--251707_at   260555_at--259090_at  
    + ... omitted several edges


La red tiene 2488 nodos y 51897 aristas. Para crear una red, en cada paso se añade un nuevo nodo. Por lo tanto, en promedio, en cada paso deberíamos añadir el número de aristas producto de la división entre el número de aristas y el número de nodos. Ya que se debe añadir un número entero de aristas por iteración, se divide sin decimales. El cociente será el número de aristas a añadir por iteración para todos los nodos excepto el último, al que se añadirán el número de aristas restantes. La función barabasi.game recibe n (numero de vertices) y el vector numérico con el número de aristas a añadir en cada paso. Para 1000 iteraciones:


```R
number.of.added.edges <- c(rep(51897%/%2488,2487), 51897 - (51897%/%2488)*2487)

clustering.coefficients <- vector(length=1000)

for(i in 1:1000)
{
  # print(i)
  random.scale.free.graph <- barabasi.game(n=2488,out.seq=number.of.added.edges,directed=FALSE)
  clustering.coefficients[i] <- transitivity(random.scale.free.graph)
}

network.clustering.coefficient <- transitivity(gene.coexpression.network.88,type = "global")

sum(clustering.coefficients > network.clustering.coefficient) / 1000
```


0


Como es 0, se rechaza que esta red haya podido ser creada de forma aleatoria, y se puede decir que la red es *libre de escala* y *de mundo pequeño*. Como tal, tendrá las siguientes propiedades: 

* Existirán nodos concentradores (hubs).
* La propagación de la información será rápida porque presentan diámetros, radios y longitudes medias de caminos mínimos pequeños en comparación con redes aleatorias.
* Son robustas a fallos o ataques aleatorios
* Son frágiles a ataques dirigidos

En la siguiente red, se muestra para cada nodo su transitividad en función de su color (azul 0, blanco 0.5 y rojo 1), demostrando que una buena parte de ellos tienen un buen número de conexiones al ser una red libre de escala:

<img src="00_Images/jta_transitivity.png" alt="figura1" height="720" width="1080">

También se puede representar el grado de los nodos, lo que sería indicativo de la conectividad de cada uno de los nodos. Para ello, se representó la red con un gradiente de color desde azul, cuando el grado de los nodos es 0 a rojo, cuando el grado de los nodos es el máximo (1572, correspondiente a 251635_at, la anotación AT3G57510 o ADPG1, una poligalacturonasa involucrada en la dihiscencia de las anteras).

<img src="00_Images/jta_degree.png" alt="figura1" height="720" width="1080">

### 3.3. Determinación de los hubs de la red

Los hubs son aquellos nodos cuya puntuación de hub está fuera del percentil 95. Se generaron listas de atributos para su posterior análisis de enriquecimiento en términos de ontología génica. Desafortunadamente, los hubs no mostraron enriquecimiento en términos de ontología génica. Para poder comprender qué función común podrían tener otros nodos concentradores que no entraran dentro del percentil 95, se tomaron los nodos que tenían una puntuación de hub inferior, tomando entonces todos los valores por encima del percentil 75:


```R
network.hub.scores <- hub.score(gene.coexpression.network.88)
hub.score.attributes <-network.hub.scores[["vector"]]
write.table(hub.score.attributes,file="hub_score_attributes_2.txt")

quatile.75 <- quantile(hub.score.attributes,prob=0.75)
hubs.values <- hub.score.attributes[hub.score.attributes > quatile.75]
hubs.names <- names(hubs.values)

length(hubs.names)

suppressMessages(library(annaffy))

hubs.table <- aafTableAnn(hubs.names, "ath1121501.db", aaf.handler())
saveHTML(hubs.table, file="hubs_table.html")
saveText(hubs.table, file="hubs_table.txt")
```


622


    Loading required package: ath1121501.db
    
    Loading required package: org.At.tair.db
    
    
    
    
    
    Warning message in chkPkgs(chip):
    “The ath1121501.db package does not appear to contain annotation data.”
    Warning message in result_fetch(res@ptr, n = n):
    “SQL statements must be issued with dbExecute() or dbSendStatement() instead of dbGetQuery() or dbSendQuery().”
    Warning message in result_fetch(res@ptr, n = n):
    “SQL statements must be issued with dbExecute() or dbSendStatement() instead of dbGetQuery() or dbSendQuery().”


El análisis de enriquecimiento de términos de ontología génica realizado con [GOrilla](http://cbl-gorilla.cs.technion.ac.il/) dio como resultado un somero enriquecimiento en genes relacionados con la respiración anaeróbica, con genes de respuesta a hipoxia (entre ellos, una galactosa oxidasa) y con procesos catabólicos de drogas, entre los que se incluyen varias peroxidasas, glucosil hidrolasas y proteínas parecidas a quitinasas (Tabla 1). Particularmente llamativa es la expresión de HRA1 (AT3G10040), un factor de transcripción atenuador de la respuesta por hipoxia. 

En los próximos apartados, se realiza una búsqueda de grupos o clústers de genes que describen la estructura subyacente de la red mediante técnicas de clustering.

<font size="2">**Tabla 1.** Resultado del enriquecimiento en términos de ontología génica de los nodos concentradores (percentil 75) en la red generada con un umbral de 0.88.</font>

<img src="00_Images/Tabla1.png" alt="figura1" height="720" width="900">

En la red siguiente, se hace analiza la puntuación de hub de cada uno de los nodos (azul 0, blanco 0.5 y rojo 1) y se hace un zoom a la región de la red que aglutina un mayor número de hubs. 

<img src="00_Images/jta_hubs_pro.png" alt="figura1" height="720" width="1080">

## 4. Análisis de clustering
<a class="anchor" id="clustering"></a>

La técnica de clustering identifica de forma automática agrupaciones o clústeres de elementos de acuerdo a una medida de similitud entre ellos. Dado un conjunto de genes y valores de similitud entre ellos (determinados con la matriz de coexpresión) se separan los genes en grupos en los que se maximiza la similitud intragrupo y se minimiza la similitud intergrupo.

En nuestro estudio, usaremos la correlación de Pearson, donde la similitud (o distancia) la define la diferencia entre 1 y la matriz de correlación (denominada gene.correlation). Según esto, un valor de distancia de 0 implica un comportamiento idéntico y un comportamiento de 2 un comportamiento simétrico. Con ello, se obtendrán diferentos grupos de genes, que serán analizados para conocer su función a nivel global y para poder determinar el efecto que la sobreexpresión de cada uno de los factores de transcripción tiene sobre estos grupos de genes.


```R
library("cluster")

similarity.matrix <- 1 - gene.correlation
```

### 4.1. Clustering jerárquico aglomerativo

La técnica de clustering jerárquico consiste en construir un dendograma que representa las relaciones de similitud entre los distintos elementos. En el clustering jerárquico aglomerativo, se comienza con tantos clústeres como individuos y consiste en ir formando (aglomerando) grupos según su similitud. Esto tiene la ventaja frente a otras técnicas de clustering de que no es necesario especificar un número de clústeres a priori, pero tiene la gran desventaja de que los errores son acumulativos, y errores en un paso de agrupamiento se propagarán al resto del dendrograma. 

La función *hclust* del paquete [*cluster*](https://cran.r-project.org/web/packages/cluster/index.html) de CRAN es el que usa para el clustering jerárquico (30). Para ello, usa la matriz de similitudes como distancias y el método de agrupamiento de enlace promedio (UPGMA) para recalcular la matriz de distancia tras cada agrupamiento. *Cutree* permite cortar el árbol generado en el clustering jerárquico a distintas alturas para producir distintos números de clústeres.


```R
hierarchical.clustering <- hclust(as.dist(similarity.matrix),method="average")

hclust.2 <- cutree(hierarchical.clustering,k=2)
hclust.3 <- cutree(hierarchical.clustering,k=3)
hclust.4 <- cutree(hierarchical.clustering,k=4)
hclust.5 <- cutree(hierarchical.clustering,k=5)
hclust.6 <- cutree(hierarchical.clustering,k=6)
hclust.7 <- cutree(hierarchical.clustering,k=7)
hclust.8 <- cutree(hierarchical.clustering,k=8)
hclust.9 <- cutree(hierarchical.clustering,k=9)
hclust.10 <- cutree(hierarchical.clustering,k=10)
hclust.11 <- cutree(hierarchical.clustering,k=11)
```

### 4.2. Clustering de partición en torno a centroides

En el clustering de partición en torno a centroides (PAM) se realiza una distribución de los elementos entre un número prefijado de clústeres o grupos. En lugar de construir un árbol como en el hclust, el objetivo en PAM consiste en agrupar los elementos en torno a elementos centrales llamados centroides (los centroides se definen como aquellos que minimizan la distancia a los puntos de su mismo clúster). En un primer paso, se seleccionan los centroides aleatoriamente y, en un segundo paso, se crean k clústeres asignando cada elemento al centroide más cercano. En el tercer paso se calcula cuál es el centroide para cada clúster. Se vuelve al segundo paso de forma iterativa hasta que bien no haya cambio en los clústeres o bien se alcance un número de iteraciones. Así, se evita la propagación de errores y permite obtener el elemento más central de cada clúster, aunque sea necesario especificar de antemano el número de clústeres a formar.

La función que realiza el clustering de partición en torno a centroides se llama *pam*, de *cluster* (1). Recibe como entrada la matriz de similitudes a usar como distancia (as.dist) y el número de clústeres a generar:


```R
pam.2 <- pam(as.dist(similarity.matrix),k=2,diss=TRUE)
pam.3 <- pam(as.dist(similarity.matrix),k=3,diss=TRUE)
pam.4 <- pam(as.dist(similarity.matrix),k=4,diss=TRUE)
pam.5 <- pam(as.dist(similarity.matrix),k=5,diss=TRUE)
pam.6 <- pam(as.dist(similarity.matrix),k=6,diss=TRUE)
pam.7 <- pam(as.dist(similarity.matrix),k=7,diss=TRUE)
pam.8 <- pam(as.dist(similarity.matrix),k=8,diss=TRUE)
pam.9 <- pam(as.dist(similarity.matrix),k=9,diss=TRUE)
pam.10 <- pam(as.dist(similarity.matrix),k=10,diss=TRUE)
pam.11 <- pam(as.dist(similarity.matrix),k=11,diss=TRUE)
```

El objetivo general de las técnicas de clustering es identificar grupos o clústeres compactos. Es decir, clusteres con una similitud intraclúster alta y una similitud interclúster baja. Esto, puede visualizarse de forma intuitiva mediante la silueta de los clústeres. La función *silhouette* nos permite calcular la silueta de un clustering, lo que sirve de medida para la bondad de dicho clustering.


```R
sil2 <- silhouette(hclust.2,dist=similarity.matrix)
sil3 <- silhouette(hclust.3,dist=similarity.matrix)
sil4 <- silhouette(hclust.4,dist=similarity.matrix)
sil5 <- silhouette(hclust.5,dist=similarity.matrix)
sil6 <- silhouette(hclust.6,dist=similarity.matrix)
sil7 <- silhouette(hclust.7,dist=similarity.matrix)
sil8 <- silhouette(hclust.8,dist=similarity.matrix)
sil9 <- silhouette(hclust.9,dist=similarity.matrix)
sil10 <- silhouette(hclust.10,dist=similarity.matrix)
sil11 <- silhouette(hclust.11,dist=similarity.matrix)

hclust.sil.values <- c(summary(sil2)[["avg.width"]],summary(sil3)[["avg.width"]],summary(sil4)[["avg.width"]],summary(sil5)[["avg.width"]],summary(sil6)[["avg.width"]],summary(sil7)[["avg.width"]],summary(sil8)[["avg.width"]],summary(sil9)[["avg.width"]],summary(sil10)[["avg.width"]],summary(sil11)[["avg.width"]])

sil2 <- silhouette(pam.2)
sil3 <- silhouette(pam.3)
sil4 <- silhouette(pam.4)
sil5 <- silhouette(pam.5)
sil6 <- silhouette(pam.6)
sil7 <- silhouette(pam.7)
sil8 <- silhouette(pam.8)
sil9 <- silhouette(pam.9)
sil10 <- silhouette(pam.10)
sil11 <- silhouette(pam.11)

pam.sil.values <- c(summary(sil2)[["avg.width"]],summary(sil3)[["avg.width"]],summary(sil4)[["avg.width"]],summary(sil5)[["avg.width"]],summary(sil6)[["avg.width"]],summary(sil7)[["avg.width"]],summary(sil8)[["avg.width"]],summary(sil9)[["avg.width"]],summary(sil10)[["avg.width"]],summary(sil11)[["avg.width"]])
```

A continuación, se representa para los dos métodos de clustering (jerárquico y PAM), y para diferentes números de clústeres la silueta correspondiente para elegir la mejor combinación de método de clustering y el número de clústeres.


```R
plot(2:11,pam.sil.values,type="o",col="blue",pch=0,ylim=c(0.2,0.5),xlab="Number of clusters",ylab="Silhouette",lwd=3)
lines(2:11,hclust.sil.values,type="o",col="red",pch=1,xlab="",ylab="",lwd=3)
legend("topright",legend=c("PAM","HCLUST"),col=c("blue","red"),pch=c(0,1),lwd=3)
```


![png](output_39_0.png)


El método PAM con 4, 5 o 6 clústeres sugiere una bondad parecida. Para comprobar cuál de los tres  es más adecuado, haremos una visualización de la silueta:


```R
sil4p <- silhouette(pam.4)
sil5p <- silhouette(pam.5)
sil6p <- silhouette(pam.6)

par(mfrow=c(1,3))
plot(sil4p,border="blue",main="PAM 4 CLUSTERS")
plot(sil5p,border="blue",main="PAM 5 CLUSTERS")
plot(sil6p,border="blue",main="PAM 6 CLUSTERS")
```


![png](output_41_0.png)


Se escogerán 4 clústeres, ya que presenta una silueta algo mejor que los 5 y 6 clústeres. Los clústeres parecen bastante compactos. En el siguiente apartado, se hace un estudio de enriquecimiento de términos de ontología de genes de cada uno de los clústeres, con el fin de determinar una función relacionada entre ellos.


```R
suppressMessages(library(annaffy))

clustering.pam.4 <- pam.4[["clustering"]]
write.table(clustering.pam.4,file="pam_4.txt")

genes.1 <- names(clustering.pam.4[clustering.pam.4==1])
genes.2 <- names(clustering.pam.4[clustering.pam.4==2])
genes.3 <- names(clustering.pam.4[clustering.pam.4==3])
genes.4 <- names(clustering.pam.4[clustering.pam.4==4])

genes.1.list <- aafTableAnn(genes.1, "ath1121501.db", aaf.handler())
saveHTML(genes.1.list, file="gclust1.html")
saveText(genes.1.list, file="gclust1.txt")

genes.2.list <- aafTableAnn(genes.2, "ath1121501.db", aaf.handler())
saveHTML(genes.2.list, file="gclust2.html")
saveText(genes.2.list, file="gclust2.txt")

genes.3.list <- aafTableAnn(genes.3, "ath1121501.db", aaf.handler())
saveHTML(genes.3.list, file="gclust3.html")
saveText(genes.3.list, file="gclust3.txt")

genes.4.list <- aafTableAnn(genes.4, "ath1121501.db", aaf.handler())
saveHTML(genes.4.list, file="gclust4.html")
saveText(genes.4.list, file="gclust4.txt")
```

Se visualizó la red en Cytoscape, manteniendo el siguiente código de colores para: cluster 1 (azul), cluster 2 (magenta), cluster 3 (verde), cluster 4 (naranja).

<img src="./00_Images/jta_cluster.png" alt="Smiley face" height="42" width="1080">


## 5. Descripción de los clústeres identificados
<a class="anchor" id="identificados"></a>

### 5.1. Introducción

En el paso anterior, se generaron tablas donde se almacenan los genes de cada clúster con su correspondiente anotación. Estos pueden ser visualizados o descargados desde la siguiente tabla:

 .| Cluster 1 | Cluster 2 | Cluster 3 | Cluster 4
-- | -- | -- | -- | --
HTML | [Ver](gclust1.html) | [Ver](gclust2.html) | [Ver](gclust3.html) | [Ver](gclust4.html)
Text | [Descargar](gclust1.txt) | [Descargar](gclust2.txt) | [Descargar](gclust3.txt) | [Descargar](gclust4.txt)

Es necesario conocer profundamente qué genes pueden estar siendo activados por estos factores de transcripción. Para ello, haremos uso de herramientas de enriquecimiento de términos de ontología de genes. La ontología es el estudio de un objeto y sus relaciones dentro de un dominio del conocimiento. En este caso, la ontología de genes (GO, Gene Ontology) describe nuestro conocimiento del dominio biológico desde tres perspectivas:

* **Procesos biológicos.** Describen eventos que ocurren de principio a fin tales como la transducción de señales.
* **Componentes celulares.** Describe las localizaciones celulares en las que el producto génico realiza una función, bien en compartimentos celulares o en complejos macromoleculares estables en los que participe.
* **Funciones moleculares.** Describen actividades que ocurren a nivel molecular tales como actividades enzimáticas (e.g. kinasas).

"El análisis de enriquecimiento de conjunto de genes (GSEA, Gene Set Enrichment Analysis) es un método computacional que determina si un conjunto de genes muestra diferencias significativas, concordantes entre dos estados biológicos". Haciendo un test exacto de Fisher puede estudiarse si existe enriquecimiento para cada GO de interés, pudiéndose determinar si en un conjunto de genes dados hay términos de GO que aparecen con significancia estadística con respecto al conjunto de genes que representan el universo total. De forma parecida, se usará [GOrilla](http://cbl-gorilla.cs.technion.ac.il/#ref) para analizar el enriquecimiento del conjunto de genes antes mencionado (6). Se comentarán (siempre que existan), los enriquecimientos para los procesos biológicos, funciones moleculares y componentes celulares.

### 5.2. Clúster 1

<img src="./00_Images/gclust1_process.png" alt="Smiley face" height="42" width="420">


En el enriquecimiento biológico, se encontró que los 11 factores de transcripción afectan a genes de síntesis de macromoléculas entre las que se encuentran muchas proteínas ribosomales e incluso proteínas de la RNA polimerasa III (AT5G09380), cuya función principal es la síntesis de RNA pequeño, así como rRNA.

<img src="./00_Images/gclust1_component.png" alt="Smiley face" height="42" width="560">

Según el enriquecimiento en componentes, se aprecia que muchos de los genes intervienen en la formación de ribosomas. De manera mucho más significativa, la expresión de estos factores de transcripción provoca un aumento de los genes que actúan en los plastos, concretamente en el cloroplasto (tanto en el estroma como en la membrana externa).

### 5.3. Clúster 2


<img src="./00_Images/gclust2_process.png" alt="Smiley face" height="42" width="175">


En el clúster 2, sólo se vio enriquecida la respuesta a estímulos, y concretamente la respuesta a químicos, como algunas proteínas de choque térmico y proteínas relacionadas con la respuesta a auxinas, que permiten el desarrollo de la planta.

### 5.4. Clúster 3

<img src="./00_Images/gclust3_component.png" alt="Smiley face" height="42" width="175">

Los transcritos del tercer cluster corresponden a proteínas extracelulares, varias relacionadas con la defensa de la planta (como defensinas como AT3G05730). También aparecen enriquecidos transcritos relacionados con la pared celular como pectinasas tales como AT4G23820. 

### 5.4. Clúster 4

<img src="./00_Images/gclust4_process.png" alt="Smiley face" height="42" width="175">

El clúster 4 es el más coexpresado de toda la red, como se ha analizado a lo largo de todo este estudio. Un estudio de su enriquecimiento en términos de ontología génica permitió conocer que estaba enriquecido el proceso catabólico de drogas, entre las que se encontraban un gran número de peroxidasas (e.g. AT2G39040, AT4G30170, AT4G26010).

### 5.5. Análisis de expresión media de los clústeres en cada condición y discusión

A continuación, para comprender en qué medida están expresados cada uno de estos clústeres cuando se sobreexpresan algunos de los 11 factores de transcripción, se representan las medias de expresión para los genes de cada clúster obtenidos a partir de la matriz traspuesta de mean.expression (diff.expr) como se explicó a lo largo del texto. Así, se podrán representar los niveles de expresión medios de cada clúster para cada genotipo.


```R
suppressMessages(library(ggplot2))

expr.cluster1.pam4 <- diff.expr[,genes.1]
expr.cluster2.pam4 <- diff.expr[,genes.2]
expr.cluster3.pam4 <- diff.expr[,genes.3]
expr.cluster4.pam4 <- diff.expr[,genes.4]

mean.profile.cluster1.pam4 <- rowMeans(expr.cluster1.pam4)
mean.profile.cluster2.pam4 <- rowMeans(expr.cluster2.pam4)
mean.profile.cluster3.pam4 <- rowMeans(expr.cluster3.pam4)
mean.profile.cluster4.pam4 <- rowMeans(expr.cluster4.pam4)

mean.vector <- c(mean.profile.cluster1.pam4, mean.profile.cluster2.pam4, mean.profile.cluster3.pam4,
                 mean.profile.cluster4.pam4)

mean.vector.for.plot <- mean.vector[order(factor(names(mean.vector), levels=tolower(c("col0", "crf2", 
                                                                                      "crf3", "czf1", "dear","erf5", "hsfc1",
                                                                                      "myb44","myb73","rav1","zat10","zf"))))]
names.means.vector <- names(mean.vector.for.plot)
cluster.for.plot <- rep(c("Clúster 1","Clúster 2","Clúster 3","Clúster 4"),4)
data <- data.frame(mean.vector,names.means.vector,cluster.for.plot)



ggplot(data, aes(fill=names.means.vector, y=mean.vector.for.plot, x=cluster.for.plot)) + 
  geom_bar(position="dodge" , stat="identity") + coord_cartesian(ylim = c(4, 8)) +
  scale_fill_brewer(palette="Paired") +labs (title = "Expresión media por clúster",
                                          subtitle = "Medida para cada genotipo", 
                                          fill = "Genotipo") +
  xlab("Clústeres") + ylab("Niveles de expresión media (a.u.)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                      panel.background = element_blank(), axis.line = element_line(colour = "black"))

```


![png](output_50_0.png)


Como se vio en el punto 2.1., y de acuerdo con los autores del artículo en revisión, el efecto sobre el transcriptoma de la mutación para sobreexpresar CRF2, DEAR, ERF5 y MYB73 no parece sustancial (1). Según este análisis, tampoco lo fue para RAV1 o ZAT10. No obstante, sí se ve un cambio sustancial para CZF1 o para HSFC1, y un efecto atenuado para CRF3, MYB44 y ZF.

En el primer clúster de genes, relacionados con procesos biosintéticos que influyen sobre todo a nivel del cloroplasto, se observó un sustancial aumento de la expresión con la sobreexpresión del factor CZF1 (*C. albicans* Zinc Finger protein), un factor de transcripción relacionado (aparte de en la respuesta al frío) con la respuesta a quitina presente en las paredes celulares de los hongos. El aumento en la síntesis de proteínas podría estar relacionado con la síntesis de genes de respuesta primaria a estos patógenos.

Dhamgaye et al. (2012) encontraron que, en *Candida albicans*, CZF1 actúa como un represor de la síntesis de $\beta$ glucanos de la pared celular, de tal forma que su ausencia causa susceptibilidad a varias drogas (Figura 2) (11). El hecho de actuar como un inhibidor de síntesis de glucanos es apreciable porque las células $\Delta$*CZF1*/$\Delta$*CZF1* son más resistentes que $\Delta$*CZF1*/CZF1, lo que indica haploinsuficiencia en estas últimas. Aplicado a nuestro estudio, y tomando la consideración de que se tratan de taxones diferentes, la sobreexpresión de CZF1 podría haber causado una represión más fuerte y sin control de los genes codificante para $\beta$ glucanos, resultando en un descenso de la resistencia disminuida a la congelación (Figura 3) (15). 

Por otro lado, se ve aprecia una disminución en los genes del clúster 4 con la sobreexpresión de CZF1. Esto podría indicar que CZF1 no solo reprimiría la expresión de los genes de síntesis de glucanos, sino también la de genes involucrado con el catabolismo de drogas. Esto es una mera interpretación de estos resultados, que podrían no reflejar la compleja realidad de la respuesta a diferentes estresores. Contrariamente a esta idea, en Tiwari et al. (2015) se describe, en respuesta a la aclimatación fría, una tolerancia aumentada a un insecticida en *Diaphorina citri* , conocido vulgarmente también como el Psílido asiático de los cítricos y considerada una de las plagas mas importante en las plantaciones de cítricos a nivel mundial (12). Estudios futuros deberán establecer una sinergia entre los mecanismos moleculares por los que CZF1 afecta a genes de síntesis de pared (11) así como a cerca de 50 genes de respuesta a frío (15) y a genes del catabolismo de drogas.

<img src="./00_Images/figura2.png" alt="Smiley face" height="42" width="600">

<span style="font-size:0.85em">**Figura 2.** Ensayos de crecimiento para probar la sensibilidad de los mutantes $\Delta$*CZF1*/$\Delta$*CZF1* y $\Delta$*CZF1*/CZF1 a drogas asociadas al SDS, FK520 (un inhibidor específico de calcineurina) y al Congo Red (que actúa en la estructura de glucanos de la pared celular). Tomado de Tiwari, S., Liu, B., Mann, R., Killiny, N. and Stelinski, L. (2015). Effects of Cold-Acclimation, Pathogen Infection, and Varying Temperatures on Insecticide Susceptibility, Feeding, and Detoxifying Enzyme Levels inDiaphorina citri (Hemiptera: Liviidae). *Florida Entomologist*, 98(3), pp.870-879.</span> 

Los descubrimientos previamente expuestos aplican a taxones muy diferentes de la planta del artículo en revisión, (*Arabidopsis thaliana*) pero se han encontrado casos de resistencia a herbicidas en plantas "voluntarias" de cultivo de colza (*Brassica napus*), que pertenecen a la misma familia (*brassicaceae*) que *Arabidopsis thaliana* (13). Los autores describen una pequeña reducción en la eficacia de los herbicidas devido a la tolerancia relacionada a la aclimatación fría o a otros factores que pueden dar lugar a la supervivencia de plantas de canola (13). Notablemente, en el proceso de aclimatación fría, la mutación de los genes CBF ($\delta$cbf123) no alteró la expresión de CZF1 (14), así como la sobreexpresión de CZF1 no alteró la expresión de los genes CBF, pero sí participó en la expresión de genes pertenecientes al regulón de CBF2; esto es, partipan en la coexpresión de genes del regulón de CBF2 sin alterar la expresión de los factores CBF (15). 

<img src="./00_Images/figura3.png" alt="Smiley face" height="42" width="300">

<span style="font-size:0.85em">**Figura 3.** Experimento de fuga de electrolitos. Estos experimentos reflejan daño a las membranas celulares, ya que son una medida de la permeabilidad de esta. Una mayor fuga de electrolitos indica una mayor permeabilidad de membrana y una reducción de la tolerancia celular al cambio de temperatura. Se muestran los resultados para la línea silvestre de control (línea continua negra) y para dos líneas transgénicas independientes que sobreexpresaban CZF1 (líneas discontinuas). Tomado de Park, S., Lee, C., Doherty, C., Gilmour, S., Kim, Y. and Thomashow, M. (2015). Regulation of the Arabidopsis CBF regulon by a complex low-temperature regulatory network. *The Plant Journal*, 82(2), pp.193-207.</span> 

CZF1 también aumenta ligeramente la expresión de proteínas del segundo clúster. No obstante, este efecto es muy ligero si se compara con el efecto que tiene la sobreexpresión del gen *hsfc1* (factor C-1 de transcripción de respuesta al estrés por calor). Este, que difiere estructuralmente y filogenéticamente (Figura 4) de otros factores de respuesta a choque térmico, se une específicamente también a elementos de choque térmico (HSE), de secuencia 5'-AGAAnnTTCT-3'(16). Los fatores de transcripción de estrés térmico son componentes finales de la transducción de múltiples señales que desencadenan la activación de genes de respuesta tanto a estrés térmico como a un gran número de estresantes químicos (16). En este contexto, es coherente que un aumento de la expresión del factor HSFC1 se correlacione con un aumento de expresión de genes relacionados con la respuesta a químicos. Más aún cuando es conocido que los factores de choque térmico juegan papeles importantes en la señalización mediada por especies reactivas de oxígeno (ROS), generadas en una gran variedad de condiciones de estrés (17).

<img src="./00_Images/philogeny.jpeg" alt="Smiley face" height="42" width="480">

<span style="font-size:0.85em">**Figura 4.** Árbol filogenético de los diferentes factores de respuesta a estrés térmico de diferentes especies de plantas. Tomado de Nover, L., Bharti, K., Döring, P., Mishra, S., Ganguli, A. and Scharf, K. (2001). Arabidopsis and the heat stress transcription factor world: how many heat stress transcription factors do we need?. Cell Stress & Chaperones, 6(3), p.177. </span> 

Con respecto a la sobreexpresión en el tercer clúster, relacionado con proteínas extracelulares, se ve que HSFC1 reduce mucho la expresión de estas. Se ha propuesto que otros factores inducibles por calor puedan afectar negativamente la transcripción de defensinas y otras proteínas extracelulares (17), pero aún queda por determinar cuáles son los mecanismos moleculares que rigen la inhibición de la expresión de estas moléculas con la sobreexpresión de HSFC1. En plantas, la congelación provoca deshidratación y estrés mecánico en la membrana plasmática, lo que puede ser uno de los factores que desencadenen la respuesta por factores de estrés como HSFC1, como se ha visto para CBF2 (18). De cualquier modo, el efecto conjunto de la sobreexpresión de HSFC1 le confirió a las plantas una mayor resistencia a las temperaturas de congelación sin aclimatación que el fenotipo silvestre. Según esto, paradógicamente, un factor de respuesta a choque térmico estaría involucrado en la aclimatación al frío. 

<img src="./00_Images/figura5.png" alt="Smiley face" height="42" width="300">

<span style="font-size:0.85em">**Figura 5.** Experimento de fuga de electrolitos. Se muestran los resultados para la línea silvestre de control (línea continua negra) y para dos líneas transgénicas independientes que sobreexpresaban CZF1 (líneas discontinuas). Tomado de Park, S., Lee, C., Doherty, C., Gilmour, S., Kim, Y. and Thomashow, M. (2015). Regulation of the Arabidopsis CBF regulon by a complex low-temperature regulatory network. *The Plant Journal*, 82(2), pp.193-207.</span> 

## 6. Conclusiones
<a class="anchor" id="conclusion"></a>

Los estudios del transcriptoma a genoma completo suponen un gran avance en biología molecular en tanto que permiten conocer los efectos de una modificación genética a nivel global y más aún cuando se combinan para conocer la convergencia a nivel transcriptómico de distintas condiciones tales como modificaciones genéticas. En esta revisión, se hizo un estudio extensivo al presentado en el artículo de Park et al. (2015) con el fin de conocer los efectos comunes de los factores de transcripción a nivel global, y no exclusivamente atendiendo a los genes de respuesta al frío. De los 11 factores de transcripción analizados (MYB73, CRF2, RAV1, CRF3, ERF5, DEAR1, MYB44, CZF1, ZAT10, ZF y HSFC1), según nuestros resultados y coherentes con los obtenidos en el artículo de referencia, no se observó un efecto sustancial en el transcriptoma para 4 (MYB73, CRF2, DEAR y ERF5). Solo se apreció un efecto atenuado para CRF3, MYB44 y ZF, y un cambio considerable para CZF1 y HSFC1. La función y los mecanismos moleculares de la transducción de señales mediadas en esta respuesta al frío parece ser muy compleja, y el análisis de la función de CZF1 y HSFC1 revela efectos muy distintos entre sí. Futuros estudios deberán estar enfocados a abarcar un mayor número de factores de transcripción con el fin de conocer todos los protagonistas de la señalización en el proceso de respuesta a temperaturas de congelación para sí reconocer los genes involucrados en la respuesta al frío. Esto, en conjunto, podría ser útil para diseñar plantas más resistentes a la congelación en vistas a paliar los efectos de las fluctuaciones de temperaturas que se esperan por el cambio climático (19).

## 7. Referencias
<a class="anchor" id="referencias"></a>


Como referencia, se tomó el artículo:

1. **Park, S., Lee, C., Doherty, C., Gilmour, S., Kim, Y. and Thomashow, M. (2015). Regulation of the Arabidopsis CBF regulon by a complex low-temperature regulatory network. *The Plant Journal*, 82(2), pp.193-207.**

Adicionalmente, se consultaron:
<font size="2">
2. Hannah, M., Heyer, A., & Hincha, D. (2005). A Global Survey of Gene Regulation during Cold Acclimation in Arabidopsis thaliana. Plos Genetics, 1(2), e26. doi: 10.1371/journal.pgen.00100262

3. Vogel, J., Zarka, D., Van Buskirk, H., Fowler, S., & Thomashow, M. (2004). Roles of the CBF2 and ZAT12 transcription factors in configuring the low temperature transcriptome of Arabidopsis. The Plant Journal, 41(2), 195-211. doi: 10.1111/j.1365-313x.2004.02288.x

4. Thomashow, M. (2010). Molecular Basis of Plant Cold Acclimation: Insights Gained from Studying the CBF Cold Response Pathway: Figure 1. Plant Physiology, 154(2), 571-577. doi: 10.1104/pp.110.161794

5. Medina, J., Bargues, M., Terol, J., Pérez-Alonso, M., & Salinas, J. (1999). The Arabidopsis CBF Gene Family Is Composed of Three Genes Encoding AP2 Domain-Containing Proteins Whose Expression Is Regulated by Low Temperature but Not by Abscisic Acid or Dehydration. Plant Physiology, 119(2), 463-470. doi: 10.1104/pp.119.2.463

6. Gilmour, S., Fowler, S., & Thomashow, M. (2004). Arabidopsis Transcriptional Activators CBF1, CBF2, and CBF3 have Matching Functional Activities. Plant Molecular Biology, 54(5), 767-781. doi: 10.1023/b:plan.0000040902.06881.d4

7. Romero-Campero, F. Introducción a la teoría de redes via *shinyapps.io*. Retrieved 2 January 2020, from https://frannetworks.shinyapps.io/red_trenes/

8. Kilian, J., Whitehead, D., Horak, J., Wanke, D., Weinl, S., & Batistic, O. et al. (2007). The AtGenExpress global stress expression data set: protocols, evaluation and model data analysis of UV-B light, drought and cold stress responses. The Plant Journal, 50(2), 347-363. doi: 10.1111/j.1365-313x.2007.03052.x

9. Podevin, N. and du Jardin, P. (2012). Possible consequences of the overlap between the CaMV 35S promoter regions in plant transformation vectors used and the viral gene VI in transgenic plants. GM Crops & Food, 3(4), pp.296-300.

10. Liu, Y., Dang, P., Liu, L., & He, C. (2019). Cold acclimation by the CBF–COR pathway in a changing climate: Lessons from Arabidopsis thaliana. Plant Cell Reports, 38(5), 511-519. doi: 10.1007/s00299-019-02376-3

11. Dhamgaye, S., Bernard, M., Lelandais, G., Sismeiro, O., Lemoine, S., Coppée, J., Le Crom, S., Prasad, R. and Devaux, F. (2012). RNA sequencing revealed novel actors of the acquisition of drug resistance in Candida albicans. BMC Genomics, 13(1), p.396.

12. Tiwari, S., Liu, B., Mann, R., Killiny, N. and Stelinski, L. (2015). Effects of Cold-Acclimation, Pathogen Infection, and Varying Temperatures on Insecticide Susceptibility, Feeding, and Detoxifying Enzyme Levels inDiaphorina citri(Hemiptera: Liviidae). Florida Entomologist, 98(3), pp.870-879.

13. Légère, A., Simard, M.-J., Johnson, E., Stevenson, F. C., Beckie, H., & Blackshaw, R. E. (2006). Control of Volunteer Canola with Herbicides: Effects of Plant Growth Stage and Cold Acclimation. Weed Technology, 20(02), 485–493. doi:10.1614/wt-05-169.1 

14. Zhao, C., Zhang, Z., Xie, S., Si, T., Li, Y., & Zhu, J. (2016). Mutational Evidence for the Critical Role of CBF Genes in Cold Acclimation in Arabidopsis. Plant Physiology, pp.00533.2016. doi: 10.1104/pp.16.00533

15. Park, S., Lee, C., Doherty, C., Gilmour, S., Kim, Y. and Thomashow, M. (2015). Regulation of the Arabidopsis CBF regulon by a complex low-temperature regulatory network. The Plant Journal, 82(2), pp.193-207.

16. Nover, L., Bharti, K., Döring, P., Mishra, S., Ganguli, A. and Scharf, K. (2001). Arabidopsis and the heat stress transcription factor world: how many heat stress transcription factors do we need?. Cell Stress & Chaperones, 6(3), p.177.

17. Kumar, M., Busch, W., Birke, H., Kemmerling, B., Nürnberger, T. and Schöffl, F. (2009). Heat Shock Factors HsfB1 and HsfB2b Are Involved in the Regulation of Pdf1.2 Expression and Pathogen Resistance in Arabidopsis. Molecular Plant, 2(1), pp.152-165.

18. Mori, K., Renhu, N., Naito, M., Nakamura, A., Shiba, H., Yamamoto, T., Suzaki, T., Iida, H. and Miura, K. (2018). Ca2+-permeable mechanosensitive channels MCA1 and MCA2 mediate cold-induced cytosolic Ca2+ increase and cold tolerance in Arabidopsis. Scientific Reports, 8(1).

19. Poutou, E., Krinner, G., Genthon, C. and de Noblet-Ducoudré, N. (2004). Role of soil freezing in future boreal climate change. Climate Dynamics, 23(6), pp.621-639.

**Transparencias del curso de biología molecular de sistemas, Francisco J. Romero Campero e Ignacio Pérez Hurtado de Mendoza, CC BY-NC-ND 3.0**

Para el análisis bioinformático, se usaron las siguientes herramientas y paquetes:

20. Carlson M (2016). ath1121501.db: Affymetrix Arabidopsis ATH1 Genome Array annotation data (chip ath1121501). R package version 3.2.3.

21. R Core Team (2019). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

22. Gautier L, Cope L, Bolstad BM, Irizarry RA (2004). “affy—analysis of Affymetrix GeneChip data at the probe level.” Bioinformatics, 20(3), 307–315. ISSN 1367-4803, doi: 10.1093/bioinformatics/btg405.

23. Miller CJ (2019). simpleaffy: Very simple high level analysis of Affymetrix data. http://www.bioconductor.org, http://bioinformatics.picr.man.ac.uk/simpleaffy/.

24. Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma powers differential expression analyses for RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7), e47. doi: 10.1093/nar/gkv007.

25. Csardi G, Nepusz T (2006). “The igraph software package for complex network research.” InterJournal, Complex Systems, 1695. http://igraph.org.

26. Shannon P, Markiel A, Ozier O, Baliga NS, Wang JT, Ramage D, Amin N, Schwikowski B, Ideker T (2003). Cytoscape: a software environment for integrated models of biomolecular interaction networks Genome Research  Nov; 13(11):2498-504

27. Eran Eden, Roy Navon, Israel Steinfeld, Doron Lipson and Zohar Yakhini (2009). GOrilla: A Tool For Discovery And Visualization of Enriched GO Terms in Ranked Gene Lists, BMC Bioinformatics , 10:48.
 
28. Eran Eden, Doron Lipson, Sivan Yogev, Zohar Yakhini (2007). Discovering Motifs in Ranked Lists of DNA sequences, PLoS Computational Biology, 3(3):e39.

29. Supek, F., Bošnjak, M., Škunca, N., & Šmuc, T. (2011). REVIGO Summarizes and Visualizes Long Lists of Gene Ontology Terms. Plos ONE, 6(7), e21800. doi: 10.1371/journal.pone.0021800

30. Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M., Hornik, K.(2019).  cluster: Cluster Analysis Basics and Extensions. R package version 2.1.0. 

</font>

**Se recoge el script completo en** **[script.r](script.r)**.
