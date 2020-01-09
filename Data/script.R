###############################################################################
#                                                                             #
#                                  TAREA 3                                    #
#                      REDES DE COEXPRESIÓN GÉNICA - SCRIPT                   #
#                                                                             #
#              Biología Molecular de Sistemas - Grado en Bioquímica           #    
#               Joaquín Tamargo Azplicueta - joatamazp@alum.us.es             #
#                                 Enero 2020                                  #
#                                                                             #
###############################################################################

################################# RESUMEN #####################################
##
## Las plantas de regiones templadas, incluyendo Arabidopsis thaliana muestran 
## un aumento en la tolerancia a la congelación en respuesta a temperaturas bajas.
## En el artículo revisado en esta tarea, los autores determinaron en primer lugar
## genes activados en respuesta a temperaturas frías de no congelación, y sobre
## ellos elaboraron una red de coexpresión génica de 11 factores de transcripción
## de primera oleada de respuesta al frío. Los resultados fueron publicados en:
##
## Park, S., Lee, C., Doherty, C., Gilmour, S., Kim, Y., & Thomashow, M. (2015). 
## Regulation of the Arabidopsis CBF regulon by a complex low-temperature regulatory
## network. The Plant Journal, 82(2), 193-207. doi: 10.1111/tpj.12796
##
## En este script se recogen las funciones usadas para la elaboración de una red
## de coexpresión génica a nivel global, y no restringida a los genes regulados
## por frío. Además, se hacen algunos comentarios al efecto. En el texto se
## amplían con detalle comentarios de las imágenes y gráficos obtenidos.

################################# ABSTRACT ####################################
##
## Plants from temperate regions (e.g. Arabidopsis thaliana), show an increase 
## in freezing tolerance in response to low nonfreezing temperatures. The authors
## of the reviewed article have firstly determined some of the genes responsible of the cold 
## acclimatation and then they determined the coexpression network of those genes, related to 11 
## transcription first-wave factors of the cold acclimatation response. 
## Their results were published on:
##
## Park, S., Lee, C., Doherty, C., Gilmour, S., Kim, Y., & Thomashow, M. (2015). 
## Regulation of the Arabidopsis CBF regulon by a complex low-temperature regulatory
## network. The Plant Journal, 82(2), 193-207. doi: 10.1111/tpj.12796
## 
## Here, I perform a comprenhensive study of the regulatory network involved 
## in the cold acclimatation not only restricted to cold-response genes, using
## the same dataset (GSE55907). This script contains all functions that were 
## neccesary for the generation and analysis of the gene coexpression network.
##
###############################################################################

########################### LECTURA DE LOS DATOS ##############################
##
## En primer lugar, se cargan los datos de GSE55907:

library(affy)

microarray.raw.data <- ReadAffy(verbose = TRUE)

########################### ANÁLISIS DE CALIDAD ###############################

## A continuación, se facilitan las instrucciones para el 
## análisis de calidad.

## library(simpleaffy)

## plot(qc(microarray.raw.data))
## boxplot(microarray.raw.data, col=rainbow(40), las=2, cex.lab=0.2)
## hist(microarray.raw.data)

## Según el análisis realizado, se aprecia un buen estado general de las
## muestras. Presentan valores de fluorescencia de las sondas semejantes 
## a pesar de aparecer como potencialmente erróneos. Los controles para 
## la degradación de muestras son negativos. De esta forma, se asume un
## buen estado de las muestras.

## A continuación, se normalizan los datos y se extraen los valores de
## expresión.

## microarray.preprocessed.data <- rma(microarray.raw.data)

## library(affyPLM)

## boxplot(microarray.preprocessed.data,col=rainbow(40), las=2, cex.lab=0.2)

expression.levels <- exprs(microarray.preprocessed.data)


##################### ANÁLISIS DE DATOS TRANSCRIPTÓMICOS ######################
############### MASIVOS Y ANÁLISIS DE EXPRESIÓN DIFERENCIAL ###################

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
head(mean.expression)

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

library(limma)

experimental.design <- model.matrix(~ -1+factor(sort(c(seq(from=1, to=12),seq(from=1, to=12)))))
colnames(experimental.design) <- experiment.name
dimnames(experimental.design)

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

length(DEG.selection.result)

######################### ANÁLISIS DE CORRELACIÓN #############################
################# ENTRE LOS PERFILES DE EXPRESIÓN GÉNICOS #####################
##
## El criterio seguido para determinar la coexpresión de dos genes se basará 
## en la correlación entre los perfiles de expresión correspondientes a los 
## distintos experimentos. Para ello, se calcula la matriz de correlaciones 
## génicas a partir de la matriz diff.expr usando la función cor.
##
## Se extraen los niveles de expresión de los genes expresados de forma diferencial
## en una matriz de expresión génica diferencial diff.expr

traspuesta<-t(mean.expression)
rownames(traspuesta) <- experiment.name
colnames(traspuesta) <- rownames(expression.levels)
diff.expr <- traspuesta[,DEG.selection.result]
dim(diff.expr)

correlation.matrix <- cor(diff.expr)

## Se obtiene una matriz cuadrada con tantas columnas y filas como genes
## diferencialmente expresados.

correlation.matrix[1:6,1:6]
dim(correlation.matrix)

######################### CONSTRUCCIÓN DE LA RED ##############################

################### Selección de la matriz de adyacencia ######################
##
## La mayoría de redes biológicas caracterizadas de forma experimental cumplen
## con la propiedad de ser redes libre de escala; es decir, que su distribución
## del grado de los nodos siga una distribución potencial negativa.
##
## Es necesario escoger un umbral de corte en la correlación para poder determinar
## los que se encuentran co-expresados de forma significativa. Para ello, se genera
## un vector que almacene unos umbrales que van a ser recorridos con un bucle
## for. Para cada valor del umbral, se obtendrá la matriz de adyacencia, se
## elaborará con ella una red y se calculará la conectividad media de esa red.
##
## Además, ya que la condición para que siga la propiedad libre de escala
## es que siga una potencial negativa, se calcularán para cada umbral 
## los valores de R2 para el ajuste lineal sobre la transformada logarítmica
## de la distribución del grado de los nodos.
##
## El umbral debe escogerse de manera que se obtengan valores de R2 altos
## y valores de conectividad media (al menos) mayores que 10.

library(igraph)

## Se crea un vector con los umbrales y vectores vacíos para almacenar las
## conectividades medias y los R2.

thresholds <- seq(from=0.70,to=0.99,by=0.01)
mean.connectivities <- vector(length=length(thresholds))
scale.free.R2 <- vector(length=length(thresholds))
gene.correlation <- correlation.matrix

par(mfrow=c(1,1))
for(i in 1:length(thresholds))
{
  print(thresholds[i])
  
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
  h <- hist(node.degrees)
  
  ## Frecuencia absoluta del grado de los nodos
  degree.frequencies <- table(node.degrees)
  
  ## Si sigue una exponencial negativa, entonces y=b*k^a
  ## De esta forma, tomando logaritmos: log(y)=a*log(k*b)
  ##
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

plot(thresholds,mean.connectivities,type="o",col="red",lwd=3,xlab="Correlation Threshold",ylab="Mean connectivity")
plot(ylim=c(0,1),thresholds,scale.free.R2,type="o",col="blue",lwd=3,xlim=c(0.70,0.99),xlab="Correlation Threshold",ylab="Scale Free Model Fit (R²)")

names(mean.connectivities) <- thresholds
names(scale.free.R2) <- thresholds

mean.connectivities
scale.free.R2


## De esta manera, se busca un compromiso entre el ajuste a 
## una red libre de escala y el mantenimiento de una conectividad 
## considerable en la red, de tal forma que se escoge 0.88, porque
## proporciona un ajuste relativamente bueno comparado con los otros
## valores a la red libre de escala (R2 = 0.74) y una alta conectividad
## (mayor de 41 vecinos):

adjacency.88 <- (gene.correlation > 0.88) & (gene.correlation < 1)
gene.coexpression.network.88 <- graph.adjacency(adjacency.88, mode="undirected")

network <- graph.adjacency(adjacency.88, mode="undirected")

network

## 2488 nodos 51897 aristas

######################## VISUALIZACIÓN DE LA RED ##############################

## La visualización de la red se realizará en Cytoscape:

write.graph(gene.coexpression.network.88,file="ath_gene_coexpression_network_88.gml",format="gml")

######################## ESTUDIO DE LA TOPOLOGÍA ##############################

#### DETERMINACIÓN DE PROPIEDAD DE RED LIBRE DE ESCALA

## 1. TEST DE KOLMOGOROV-SMIRNOV
##
## Como se señaló previamente, se estudiará si la distribución del grado de los
## nodos sigue una exponencial negativa. Para ello, con el test de 
## Kolmogorov-Smirnov se hace un contraste de hipótesis donde: la hipótesis nula
## (H0): la distribución del grado de los nodos sigue una potencial;
## la hipótesis alternativa (H1) dice que la distribución de los grados de los
## nodos no sigue una distribución potencial.

## Para este estudio se extrae la distribución del grado de los nodos y se 
## hace la comprobación de la presencia de la ley de potencia:

network.degree.distribution <- degree.distribution(gene.coexpression.network.88)
fit.scale.free <- power.law.fit(network.degree.distribution)
fit.scale.free[["KS.p"]]

## Como 0,9487 >>> alfa, no se puede descartar la hipótesis nula: la distribución
## del grado de los nodos sigue una red libre de escala.

## 2. REGRESIÓN LINEAL
##
## Se siguen los pasos utilizados anteriormente para el apartado
## de selección de la matriz de adyacencia.

# Cálculo y representación del grado de los nodos

par(mfrow=c(1,1))
network.degrees <- degree(gene.coexpression.network.88)
degree.histogram <- hist(network.degrees,
                         freq=FALSE, col="light blue",
                         xlab="Node Degrees",
                         ylab="Probability", main="Distribución del grado de los nodos")

# Calculo de la frecuencia absoluta del grado de los nodos
degree.frequencies <- table(network.degrees)

# Eliminamos nodos con grado 0 para poder aplicar log10
degree.frequencies.no.0 <- degree.frequencies[-1]

# Transformación logarítmica
log10.degrees.frequencies <- log10(degree.frequencies.no.0)
log10.node.degrees <- log10(as.numeric(names(degree.frequencies.no.0)))

# Regresión lineal

lm.r <- lm(log10.degrees.frequencies ~ log10.node.degrees)
plot(x = log10.node.degrees, y = log10.degrees.frequencies, pch=20, col = "grey", ylim=c(0,3),
     main="Regresión de la frecuencia de los grados de los nodos 
     con transformación logarítmica", xlab="log(Grado de los nodos)", ylab="log(Frecuencia de grados)")
abline(lm.r, col="blue")
summary(lm.r)

## Como p-valor <<< 0,001, hay una correlación lineal y por 
## tanto sigue una distribución potencial negativa (a=-1,01). 

#### DETERMINACIÓN DE MUNDO PEQUEÑO
##
## La condición que debe cumplir una red de mundo pequeño es que sus nodos,
## (las entidades del sistema) aunque sean muy abundantes, estén muy
## conectados y sea posible hacer caminos cortos para llegar de un
## nodo a otro cualquiera de la red. Los valores de la transitividad
## (o coeficiente de agrupamiento) y la longitud media del camino mínimo
## son principalmente los indicadores que permitirán hacer la determinación
## de que la red sea de mundo pequeño.
##
## La transitividad representa el grado de agrupamiento en torno a un
## nodo, y toma valores entre 0 (todos los vecinos del nodo de
## interés no son vecinos entre sí) y 1 (todos los nodos vecinos del 
## nodo de interés son vecinos entre sí).

network.clustering.coefficient <- transitivity(gene.coexpression.network.88,type="global")
average.path.length(gene.coexpression.network.88, directed=FALSE)

## Tiene un coeficiente de agrupamiento de 0.72 y la longitud de 
## camino mínimo media es de 6.46. 

## Para comprobar que la generación de la red no ha tenido lugar de
## forma puramente aleatoria, se estimará la probabilidad de obtener
## una red similar a la estudiada mediante el método de Albert-Lászlo 
## Barabasi, por el cual se parte de una pequeña red donde todos los
## nodos están unidos entre sí y se van conectando nuevos nodos a los 
## ya existentes, pero no de forma arbitraria, sino con una preferencia
## mayor por los nodos de alto grado. Para ello, se usará la función
## barabasi.game, que recibe como entrada el tamaño de la red a crear
## (el número de nodos) y un vector numérico que da el número de aristas
## a añadir en cada paso, y devuelve una red libre de escala.

## En primer lugar, se evalúa la red de coexpresión génica generada
## para comprobar el número de nodos y aristas de la red.

gene.coexpression.network.88

## La red tiene 2488 nodos y 51897 aristas

## En cada paso se añade un nuevo nodo. Por lo tanto, en promedio, en cada 
## paso deberíamos añadir el siguiente número de aristas: 
## numero_de_aristas/numero_de_nodos. Ya que se debe añadir un número entero 
## de aristas por iteración, se divide sin decimales. El cociente será el número
## de aristas a añadir por iteración para todos los nodos excepto el último, al
## que se añadirán el número de aristas restantes.

## Se calcula el cociente sin decimales:

51897%/%2488

## En primer lugar, se añaden 20 aristas a cada uno de los 2488 nodos. Luego, 
## habrá que añadir sobre un único nodo el resto:

51897 - (51897%/%2488)*2487

## De esta manera, sumando ambos se tienen que obtener 51897 aristas:

sum(c(rep(51897%/%2488,2487), 51897 - (51897%/%2488)*2487))

## Así, se generan grafos aleatorios y calculamos su coeficiente de agrupamiento medio

number.of.added.edges <- c(rep(51897%/%2488,2487), 51897 - (51897%/%2488)*2487)

## La función barabasi.game recibe n (numero de vertices), un vector numerico con el numero de
## aristas a añadir en cada paso. Su primer elemento es ignorado porque no se añaden aristas
## en el primer paso.

clustering.coefficients <- vector(length=1000)

for(i in 1:1000)
{
  print(i)
  random.scale.free.graph <- barabasi.game(n=2488,out.seq=number.of.added.edges,directed=FALSE)
  clustering.coefficients[i] <- transitivity(random.scale.free.graph)
}

network.clustering.coefficient <- transitivity(gene.coexpression.network.88,type = "global")

sum(clustering.coefficients > network.clustering.coefficient) / 1000

## Como es 0, se puede rechazar que la red haya sido generada de forma
## aleatoria, y concluimos que la red es de mundo pequeño.


################### GENERACIÓN DE TABLAS PARA CYTOSCAPE #######################

#### TRANSITIVIDAD (Coeficiente de agrupamiento)

node.transitivity <- transitivity(gene.coexpression.network.88, type = "local", isolates="zero")
names(node.transitivity) <- names(V(gene.coexpression.network.88))
write.table(node.transitivity, file="node_transitivity.txt")

#### IDENTIFICACIÓN DE NODOS CONCENTRADORES

## En las redes libres de escala, la probabilidad de encontrar un nodo con un 
## bajo grado de forma aleatoria es alto, mientras que la probabilidad de encontrar
## un nodo con un alto grado de forma aleatoria es bajo pero no despreciable.
## En estas redes, estos nodos concentradores con muchas aristas se conocen
## como hubs, y la eliminación de estos tiene un efecto muy fuerte sobre la 
## red.

network.hub.scores <- hub.score(gene.coexpression.network.88)
hub.score.attributes <-network.hub.scores[["vector"]]
write.table(hub.score.attributes,file="hub_score_attributes_2.txt")

#### GRADO DE LOS NODOS

network.degree.distribution <- degree(gene.coexpression.network.88)
write.table(network.degree.distribution, file="network_degree.txt")

######################## IDENTIFICACIÓN DE LOS HUBS ###########################

## Los HUBS los nodos cuya puntuación de hub está fuera del percentil 95. 
## A continuación se generan listas de atributos para su posterior análisis
## de enriquecimiento en términos de ontologíagénica.

quatile.95 <- quantile(hub.score.attributes,prob=0.95)
hubs.values <- hub.score.attributes[hub.score.attributes > quatile.95]
hubs.names <- names(hubs.values)
length(hubs.names)

## Desafortunadamente, los HUBS no mostraron enriquecimiento en términos de 
## ontología génica. Para poder comprender qué función común podrían 
## tener otros nodos concentradores que no entraran dentro del percentil 
## 95, se tomaron los nodos que tenían una puntuación de hub inferior,
## tomando entonces todos los valores por encima del percentil 75.

quatile.75 <- quantile(hub.score.attributes,prob=0.75)
hubs.values <- hub.score.attributes[hub.score.attributes > quatile.75]
hubs.names <- names(hubs.values)

length(hubs.names)

library(annaffy)

hubs.table <- aafTableAnn(hubs.names, "ath1121501.db", aaf.handler())
saveHTML(hubs.table, file="hubs_table.html")
saveText(hubs.table, file="hubs_table.txt")


############################## CLUSTERING ####################################

## La técnica de clustering identifica de forma automática agrupaciones o 
## clústeres de elementos de acuerdo a una medida de similitud entre ellos.
## Dado un conjunto de genes y valores de similitud entre ellos (coexpresión)
## se determinan grupos en los que se maximiza la similitud intragrupo y se 
## minimiza la similitud intergrupo.

## Todas las técnicas de clustering buscan minimizar la distancia intracluster y 
## maximizar las distancias intercluster. En nuestro estudio, usaremos el la
## correlación de Pearson, donde la distancia la define la diferencia entre
## 1 y la matriz de correlación (que recibe el nombre de gene.correlation).
## Según esto, un valor de distancia de 0 implica un comportamiento idéntico
## y un comportamiento de 2 un comportamiento simétrico.

library("cluster")

similarity.matrix <- 1 - gene.correlation


#### CLUSTERING JERÁRQUICO

## La función hclust usa la matriz de similitudes como distancias y el método promedio
## para recalcular distancias calcula el clustering jerárquico.

hierarchical.clustering <- hclust(as.dist(similarity.matrix),method="average")

## La función cutree permite cortar el árbol generado en el clustering jerárquico a distintas
## alturas para producir distintos números de clústeres.

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

#### CLUSTERING DE PARTICIÓN EN TORNO A CENTROIDES

## En el clustering de partición en torno a centroides (PAM) se realiza
## una distribución de los elementos entre un número prefijado de
## clústeres o grupos. Esta técnica recibe como dato de entrada 
## el número de clústers a formar (k) además de los elementos a 
## clasificar y la matriz de similitudes como distancias.
## 
## En lugar de construir un árbol como en el hclust, el objetivo en
## PAM consiste en agrupar los elementos en torno a elementos 
## centrales llamados centroides.

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

## La función silhouette nos permite calcular la silueta de un clustering que sirve de medida
## para la bondad de dicho clustering.

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

plot(sil5,border="blue", main = "SILUETA PARA PAM5")

pam.sil.values <- c(summary(sil2)[["avg.width"]],summary(sil3)[["avg.width"]],summary(sil4)[["avg.width"]],summary(sil5)[["avg.width"]],summary(sil6)[["avg.width"]],summary(sil7)[["avg.width"]],summary(sil8)[["avg.width"]],summary(sil9)[["avg.width"]],summary(sil10)[["avg.width"]],summary(sil11)[["avg.width"]])

## Se representa para los dos métodos de clustering (jerárquico y PAM), y para diferentes números
## de clústeres la silueta correspondiente para elegir la mejor combinación de método de 
## clustering y el número de clústeres.

plot(2:11,pam.sil.values,type="o",col="blue",pch=0,ylim=c(0,0.8),xlab="Number of clusters",ylab="Silhouette",lwd=3)
lines(2:11,hclust.sil.values,type="o",col="red",pch=1,xlab="",ylab="",lwd=3)
legend("topright",legend=c("PAM","HCLUST"),col=c("blue","red"),pch=c(0,1),lwd=3)

## Se escogerán 4, 5 o 6 clústeres. Para comprobar cuál de los tres  es más adecuado, haremos una
## visualización de la silueta:

sil4p <- silhouette(pam.4)
sil5p <- silhouette(pam.5)
sil6p <- silhouette(pam.6)

par(mfrow=c(1,3))

plot(sil4p,border="blue",main="PAM 4 CLUSTERS")
plot(sil5p,border="blue",main="PAM 5 CLUSTERS")
plot(sil6p,border="blue",main="PAM 6 CLUSTERS")

## La silueta tiene mejor aspecto para el método PAM4. A continuación
## se guardan los datos en una tabla que puede introducirse en Cytoscape.

clustering.pam.4 <- pam.4[["clustering"]]
write.table(clustering.pam.4,file="pam_4.txt")

## Se almacenan los genes de cada cluster en un HTML y un txt

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


## Se obtienen las medias de expresión para los genes de cada clúster a partir
## de la matriz traspuesta de mean.expression (diff.expr). Así, se podrán 
## representar los niveles de expresión medios de cada clúster para cada genotipo.

dim(diff.expr)

expr.cluster1.pam4 <- diff.expr[,genes.1]
expr.cluster2.pam4 <- diff.expr[,genes.2]
expr.cluster3.pam4 <- diff.expr[,genes.3]
expr.cluster4.pam4 <- diff.expr[,genes.4]

mean.profile.cluster1.pam4 <- rowMeans(expr.cluster1.pam4)
mean.profile.cluster2.pam4 <- rowMeans(expr.cluster2.pam4)
mean.profile.cluster3.pam4 <- rowMeans(expr.cluster3.pam4)
mean.profile.cluster4.pam4 <- rowMeans(expr.cluster4.pam4)

###

library(ggplot2)

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


############################## REFERENCIAS ####################################
##
##  1. Park, S., Lee, C., Doherty, C., Gilmour, S., Kim, Y. and Thomashow, M. (2015).
##     Regulation of the Arabidopsis CBF regulon by a complex low-temperature regulatory 
##     network. The Plant Journal*, 82(2), pp.193-207.**
##  
##  2. Romero-Campero, F. Introducción a la teoría de redes via *shinyapps.io*. 
##     Retrieved 2 January 2020, from https://frannetworks.shinyapps.io/red_trenes/
##
##  3. Transparencias del curso de biología molecular de sistemas, Francisco J. Romero
##     Campero e Ignacio Pérez Hurtado de Mendoza, CC BY-NC-ND 3.0**
##
##  4. Carlson M (2016). ath1121501.db: Affymetrix Arabidopsis ATH1 Genome Array annotation 
##     data (chip ath1121501). R package version 3.2.3.
##
##  5. R Core Team (2019). R: A language and environment for statistical computing. R 
##     Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
##
##  6. Gautier L, Cope L, Bolstad BM, Irizarry RA (2004). “affy—analysis of Affymetrix 
##     GeneChip data at the probe level.” Bioinformatics, 20(3), 307–315. ISSN 1367-4803, 
##     doi: 10.1093/bioinformatics/btg405.
##
##  7. Miller CJ (2019). simpleaffy: Very simple high level analysis of Affymetrix data.
##     http://www.bioconductor.org, http://bioinformatics.picr.man.ac.uk/simpleaffy/.
##
##  8. Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma powers
##     differential expression analyses for RNA-sequencing and microarray studies.” Nucleic 
##     Acids Research, 43(7), e47. doi: 10.1093/nar/gkv007.
##
##  9. Csardi G, Nepusz T (2006). “The igraph software package for complex network research.”
##     InterJournal, Complex Systems, 1695. http://igraph.org.
##
## 10. Shannon P, Markiel A, Ozier O, Baliga NS, Wang JT, Ramage D, Amin N, Schwikowski B,
##     Ideker T. Cytoscape: a software environment for integrated models of biomolecular 
##     interaction networks Genome Research 2003 Nov; 13(11):2498-504
##
## 11. Eran Eden, Roy Navon, Israel Steinfeld, Doron Lipson and Zohar Yakhini. "GOrilla:
##     A Tool For Discovery And Visualization of Enriched GO Terms in Ranked Gene Lists", 
##     BMC Bioinformatics 2009, 10:48.
##
## 12. Eran Eden, Doron Lipson, Sivan Yogev, Zohar Yakhini. "Discovering Motifs in Ranked 
##     Lists of DNA sequences", PLoS Computational Biology, 3(3):e39, 2007.
## 13. Supek, F., Bošnjak, M., Škunca, N., & Šmuc, T. (2011). REVIGO Summarizes and Visualizes
##     Long Lists of Gene Ontology Terms. Plos ONE, 6(7), e21800. doi: 10.1371/journal.pone.0021800
##
#######################################################################################################
