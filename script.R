


#                         convegno AIE di Autunno 2021

#                                 Workshop 
# "Introduzione all'approccio esposomico in epidemiologia con il software R"

#                     Chiara Moccia & Antonio D'Errico  
# Dottorandi Dip.to Scienze Mediche, Epidemiologia dei Tumori, Università di Torino
#                     

#Install package rexposome>>>>eseguire i due comandi successivi solo alla prima esecuzione dello script 
install.packages("devtools")
devtools::install_github("isglobal-brge/rexposome")
# title = {rexposome: Exposome exploration and outcome data analysis},
# author = {Carles Hernandez-Ferrer and Juan R. Gonzalez and Xavier Escribà-Montagut},
# year = {2021},
# note = {R package version 1.15.0},

###Enter one or more numbers, or an empty line to skip updates:  ---> premere invio per non installare aggiornamenti

rm(list=ls(all=TRUE))
library(rexposome)
library(ggplot2)


## inserire il percorso alla cartella
path <- file.path("F:/EXPOSOME-workshop-AIE2021--main")

# rexposome richiede in input tre diversi dataset:  
# 1.	Description Data
# 2.	Exposure Data
# 3.	Phenotype Data

description <- file.path(path, "descriptions.csv")
phenotype <- file.path(path, "phenotypes.csv")
exposures <- file.path(path, "exposures.csv")

dd <- read.csv(description, header = TRUE)
View(dd)
ee <- read.csv(exposures, header = TRUE)
View(ee)
pp <- read.csv(phenotype, header = TRUE)
View(pp)

# Il description data file avrà sulle righe le esposizioni

## ----assegniamo nome delle righe a description file---------------------------
rownames(dd) <- dd$Exposure
dd$Exposure <- NULL


# L'exposures data file avrà una colonna per ogni esposizione definita nel description file
# e sulle righe gli individui

## ----assegniamo nome delle righe a exposures file-----------------------------
rownames(ee) <- ee$idnum
ee$idnum <- NULL



# Il phenotype data file avrà una colonna per ogni fenotipo/outcome
# e sulle righe gli stessi individui dell'exposures data file

## ----assegniamo nome delle righe a phenotype file-----------------------------
rownames(pp) <- pp$idnum
pp$idnum <- NULL
View(pp)
## ----creiamo l'oggetto della classe ExposomeSet, con i 3 dataframe caricati sopra
exp <- loadExposome(
  exposures = ee, 
  description = dd, 
  phenotype = pp,
  description.famCol = "Family" #indica la colonna contenente la famiglia di esposizioni nel description file
)
exp   

### Funzioni utili per accedere alle info contenute nell'oggetto exp 

## ----identificativi dei soggetti----------------------------------------------
head(sampleNames(exp))

## ----nomi delle esposizioni---------------------------------------------------
head(exposureNames(exp))

## ----nomi delle famiglie di esposizioni---------------------------------------
familyNames(exp)

## ----nomi dei fenotipi--------------------------------------------------------
phenotypeNames(exp)

## ----funzioni per accedere ai valori delle variabili--------------------------
# 
head(fData(exp), n = 3) #feature data (description data) 

## ----phenotype data-----------------------------------------------------------
head(pData(exp), n = 3) #descrizione dei fenotipi

## ----exposures data-----------------------------------------------------------
expos(exp)[1:10, c("Cotinine", "PM10V", "PM25", "X5cxMEPP")]






#--------------------------------Exposome Pre-processing------------------------
#                       Missing Data in Exposures and Phenotypes

## ----tabulare i missing data--------------------------------------------------
# tableMissings(exp, set = "exposures", output = "n") # p percentuali, n numero di missing
# tableMissings(exp, set = "phenotypes", output = "p")

## ----plot dei missing nelle esposizioni---------------------------------------
win.graph(10, 10, 12)
plotMissings(exp, set = "exposures")
## ----plot dei missing nei fenotipi--------------------------------------------
win.graph(10, 10, 12)
plotMissings(exp, set = "phenotypes")

#                             Missing Data imputation
# L'imputazione è effettuata con la funzione imputation

# exp_imp=imputation(exp)

#                               Exposures Normality

# ## ----Test di Shapiro-Wilk per la verifica della normalità---------------------
 nm <- normalityTest(exp,th=0.05) 
 nm
# table(nm$normality)
# 
# ## ----lista delle variabili che non hanno distribuzione normale----------------
# nm$exposure[!nm$normality]
# 
# #  Un'ispezione visiva permette di identificare la trasformazione che meglio 
# #  approssimi la normale grazie all'argomento show.trans=TRUE
 win.graph(10, 10, 12)
 plotHistogram(exp, select = "Dens", show.trans = TRUE, density = TRUE)#sqrt

 

#                               Exposures transformation
## Trasformo le variabili non normali tramite trasformazione più appropriata (sqrt, log, exp, inv..) 

exp_trans=trans(exp, sqrt, select="Dens")
exp_trans


#                              Exposures standardization
# La standardizzazione è effettuata con la funzione standardize

exp_std <- standardize(exp_trans, select ="Dens" )
exp_std




#                             Exposures Characterization 


familyNames(exp)
win.graph(10, 10, 12)
plotFamily(exp, family = "all")
win.graph(10, 10, 12)
plotFamily(exp, family = "Phthalates")

#plotFamily(exp, family = "Air Pollutants", group = "sex")
win.graph(10, 10, 12)
plotFamily(exp, family = "Air Pollutants", group = "sep_bin")
#plotFamily(exp, family = "Built Environment", group = "sep_bin")
#plotFamily(exp, family = "Metals", group = "sep_bin")
win.graph(10, 10, 12)
plotFamily(exp, family = "Cotinine", group = "sep_bin")

#plotFamily(exp, family = "PBDEs", group = "sex")
#plotFamily(exp, family = "Organochlorines", group = "sep_bin")
#plotFamily(exp, family = "Bisphenol A", group = "sex")
#plotFamily(exp, family = "Water Pollutants", group = "sep_bin")
#plotFamily(exp, family = "Phthalates", group = "sep_bin")
#plotFamily(exp, family = "Noise", group = "sep_bin")
#plotFamily(exp, family = "PFOAs", group = "sep_bin")
#plotFamily(exp, family = "Temperature", group = "sep_bin")






#                              Exposures correlation
# Analisi della correlazione tra le esposizioni, in termini di 
# esposizioni intra-famiglia o inter-famiglia.
exp_cr <- correlation(exp, use = "pairwise.complete.obs", method.cor = "pearson",warnings=FALSE)
extract(exp_cr)[1:4, 1:4]

## ---- circos plot: utile per inter-family correlation-------------------------
windows(rescale="fit")
plotCorrelation(exp_cr, type = "circos")
#Green inversamente correlato air pollutants
#Noise direttamente correlato con air pollutants

## ----correlation matrix: utile per intra-family correlation-------------------
windows(rescale="fit")
plotCorrelation(exp_cr, type = "matrix")




#                               Exposures PCA
# Una volta che le esposizioni sono state standardizzate, si può procedere con 
# l'analisi delle componenti principali usando il metodo pca. 
# 

# FAMD: pca per mixed data type
# exp_pca <- pca(exp_std) 

# PCA classica su variabili continue
exp_pca <- pca(exp,pca = T)

## ----plot_pca-----------------------------------------------------------------
win.graph(10, 10, 12)
plotPCA(exp_pca, set = "all",  show.exposures = T,
        show.samples = T)
# plotPCA(exp_pca, set = "samples",  show.exposures = T,
#         show.samples = T)
# plotPCA(exp_pca, set = "exposures",  show.exposures = T,
#         show.samples = T)

## ----plot_pca_samples---------------------------------------------------------
# plotPCA(exp_pca, set = "samples", phenotype = "sex")
# 
# plot3PCA(exp_pca, cmpX=1, cmpY=2, cmpZ=3, phenotype = "sep_bin")

win.graph(10, 10, 10)
plotEXP(exp_pca)+theme(axis.text.y = element_text(size=6.5))+ylab("") 



##                           Individuals Clustering
# L'analisi dei cluster dell'esposoma permette di raggruppare e caratterizzare 
# i soggetti sulla base dei loro profili di esposizione. 

args(clustering)

## ----hclust_function----------------------------------------------------------
hclust_data <- function(data, ...) {
  #hclust(d = dist(x = data,method = "manhattan"), method="single")
  hclust(d = dist(x = data))# "euclidean", "complete" by default
  #method in dist: the distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
  #method in hclust: the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
}

## ----hclus_k3-----------------------------------------------------------------
hclust_k3 <- function(result) {
  cutree(result, k = 3)
}


## ----expo_clustering----------------------------------------------------------
exp_c <- clustering(exp, method = hclust_data, cmethod = hclust_k3)
exp_c
exp_c@model$model$method
exp_c@model$model$dist.method


# numero individui in ognuno dei 3 cluster
table(classification(exp_c))

## ----Per ogni cluster è plottato il profilo di esposizione-----------------------------
windows(30,30,18)
plotClassification(exp_c)












#                                   EXPOSURES ASSOCIATIONS

#                 Analisi univariata: Exposome-wide Association Study (ExWAS)

# Il metodo exwas performa un exposome wide association study
# testando l'associazione tra una esposizione e un fenotipo/outcome nell'ExposomeSet.
# Attraverso il parametro formula possiamo indicare il fenotipo che vogliamo testare
# come outcome e le covariate che vogliamo includere nel modello.

# La correzione per confronti multipli nell'exwas è effettuata calcolando 
# il numero di test effettivi e la soglia di significatività (Li and Ju, https://doi.org/10.1038/sj.hdy.6800717).

bl_ew <- exwas(exp, formula = blood_pre~sex+age, family = "gaussian",tef = TRUE)
head(extract(bl_ew))

we_ew <- exwas(exp, formula = wheezing~sex+age, family = "binomial",tef = TRUE)
head(extract(we_ew))

#plotExwas: plotta per ogni esposizione il pvalue dell'associazione in un manhattan plot
clr <- rainbow(length(familyNames(exp)))
names(clr) <- familyNames(exp)
win.graph(10, 10, 12)
plotExwas(bl_ew, we_ew, color = clr,show.effective=TRUE) + ggtitle("Exposome Association Study - Univariate Approach")

#plotEffect: plotta la stima dell'effetto di ogni esposizione con il relativo intervallo di confidenza 
win.graph(10, 10, 12)
plotEffect(bl_ew)
win.graph(10, 10, 12)
plotEffect(we_ew)







#                     Associazione tra appartenenza al cluster e un fenotipo
win.graph(10, 10, 12)
plotClassification(exp_c)

cluster=as.factor(exp_c@phenoData@data$cluster)
bl_press=exp_c@phenoData@data$blood_pre
age=exp_c@phenoData@data$age
sex=exp_c@phenoData@data$sex
fit = glm(formula=bl_press ~ cluster+age+sex,family=gaussian)
summary(fit)







#                       Associazione tra Componenti Principali e un fenotipo

# Dai risultati della PCA sull'ExposomeSet, otteniamo due utili misure:
# Correlazione tra esposizioni e componenti principali
# Associazione tra fenotipi e componenti principali
win.graph(10, 10, 10)
plotEXP(exp_pca)+theme(axis.text.y = element_text(size=6.5))+ylab("") 

# PC1 altamente correlata con PM
# PC2 correlata con BDEs
# PC5 negativamente correlata con PFOAs e PCBs
# Utile per caratterizzare le PC in termini di esposizioni

# Il metodo plotPHE testa l'associazione tra i fenotipi e le PC
win.graph(10, 10, 10)
plotPHE(exp_pca)
# blood pressure associata a PC1
# sep associato a PC4
# sex associato a PC10
# age associata a PC10








##            Analisi multivariata:  Multivariate-Exposome-wide Association Study (ExWAS)
# La funzione mexwas permette di testare un'analisi di associazione multivariata
# tra l'esposoma e un outcome usando Elastic-Net regularized generalized linear models
# L'obiettivo è trovare il gruppo di esposizioni maggiormente associate con l'outcome.

bl_mew <- mexwas(exp, phenotype = "blood_pre", family = "gaussian")
we_mew <- mexwas(exp, phenotype = "wheezing", family = "binomial")

# Con il metodo plotExwas possiamo plottare il coefficiente per ogni esposizione.
# Il metodo fornisce una mappa colore che evidenzia i coefficienti più interessanti.
# Le due colonne della heat map corrispondono al lambda con minimo errore medio in crossvalidazione 
# e al lambda che fornisce il modello più regolarizzato il modello più regolarizzato 
win.graph(10, 10, 12)
plotExwas(bl_mew, we_mew) + ylab("") +
  ggtitle("Exposome Association Study - Multivariate Approach")
