#################################
# Data Analysis and Exploration #
# A.A. 2019-2020                #
# Ambrosi Andrea                #
#################################

library("GEOquery")
gse  <- getGEO(GEO = 'GSE48000', filename = "/home/andrea/Scrivania/data_analysis/GSE48000_series_matrix.txt.gz") # processed data

# The actual expression data are accessible in the "exprs" ... an
# Expression Set, the generic data class that BioConductoruses
# for expression data
head(exprs(gse))
dim( exprs(gse)) # 47304 rows times 132 columns

# exprs(gse) is a matrix that we can assign to its own variable, to
# conveniently access the data rows and columns
ex <- exprs(gse)
#table(is.na(ex)) #there's no NA
dim(ex)
head(ex)
colnames(ex)

#    Analyze value distributions
#boxplot(ex) #e un casino da vedere, troppe colonne e troppi dati


#------------------------------------- w2 → w3

# PCA
pca <- prcomp(t(ex))

summary(pca)
screeplot(pca)

# draw PCA plot
grpcol <- c(rep("blue", 25), rep("red", 40), rep("green", 34), rep("yellow", 33))
plot(pca$x[ , 1], pca$x[ , 2], xlab = "PCA1", ylab = "PCA2", 
	 main = "PCA for components 1&2", type = "p", pch = 10, col = grpcol)
text(pca$x[ , 1], pca$x[ , 2], rownames(pca$x), cex = 0.75)

#------------------------------------- w3 → w4


#R CODE –K-MEANS EXAMPLE
#BiocManager::install("useful")
library("useful")

k <- 4

kmeans_result <- kmeans(t(ex), k)

table(kmeans_result$cluster)

#cluster visualization reduction using PCA 

plot(kmeans_result, data = t(ex))


#R CODE –HIERARCHICAL CLUSTERING EXAMPLE

dist_matrix <- dist(t(ex))
hc_result   <- ?hclust(dist_matrix, method = "ave")
# the agglomeration method to be used. This should be 
#(an unambiguous abbreviation of) one of "ward.D", 
#"ward.D2", "single", "complete", "average" (= UPGMA), 
#"mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).

k <- 4
groups <- cutree(hc_result, k = k)

table(groups)

plot(hc_result, hang <- -1, labels = groups)
rect.hclust(hc_result, k = 4, which = NULL, x = NULL, h = NULL, border = 2, cluster = NULL) # red boxes to show groups


#------------------------------------- w4 → w5

#CODE EXAMPLE
install.packages("BiocManager")
# load data 
BiocManager::install("ALL"); library("ALL"); data(ALL)

# keep only 30 arrays JUST for computational convenience
e.mat <- 2^(exprs(ALL)[ , c(81:110)]) 

# Genefilterpackage is very useful to preprocess data# here we remove genes whose value is not > 100 in at least 20% of the samples# filter genes on raw scale, then return to log scale; 
BiocManager::install("genefilter"); library("genefilter"); # cernita sulle righe della matrice. Non necessario sempre ma molto utile per generare e definire funzioni per selezionare le righe che vanno bene
#ffun       <- filterfun(pOverA(0.20, 100)) #elimina tutte le righe dove ...
#t.fil      <- genefilter(e.mat, ffun)
#small.eset <- log2(e.mat[t.fil, ])

small.eset <- log2(na.omit(e.mat))
dim(small.eset) # 12625 genes, 30 arrays (15 B and 15 T)
group      <- c(rep('B', 15), rep('T', 15)) # classification, in order 
#group sono le etichette dei dati da 81 a 110 scelti appositamente
# Build RF

BiocManager::install(" randomForest"); library(randomForest) 
set.seed(1234) 
print(date()) 
rf <- randomForest(x = t(small.eset), y = as.factor(group), ntree = 1000) 
print(date()) # it takes about 20 seconds

# a trivial test
predict(rf, t(small.eset[ , 1:5]))


#FACTORS
data  <- c(1, 2, 2, 3, 1, 2, 3, 3, 1, 2, 3, 3, 1) 
fdata <- factor(data) 

# >fdata[1] 1 2 2 3 1 2 3 3 1 2 3 3 1  
# Levels:1 2 3 

rdata <- factor(data, labels = c("I", "II", "III")) 

# >rdata[1] I II II III I II III III I II III III I 
# Levels:I II III

#da fare solo se interessati
#DRAWING A HEATMAP
# Look at variable importance
imp.temp <- abs(rf$importance[ , ])
t        <- order(imp.temp, decreasing = "TRUE")
plot(c(1:nrow(small.eset)), imp.temp[t], log = 'x', cex.main = 1.5     , 
	 xlab = 'gene rank'   , ylab = 'variable importance', cex.lab = 1.5,
	 pch  = 16, main = 'ALL subset results')

# plot utile nel report menzionato a lezione

# Get subset of expression values for 25 most 'important' genes
gn.imp   <- names(imp.temp)[t]
gn.25    <- gn.imp[1:25]    # vector of top 25 genes, in order
t        <- is.element(rownames(small.eset), gn.25)
sig.eset <- small.eset[t, ] # matrix of expression values, not necessarily in order

## Make a heatmap, with group differences obvious on plot
library(RColorBrewer)
hmcol              <- colorRampPalette(brewer.pal(11, "PuOr"))(256)
colnames(sig.eset) <- group # This will label the heatmapcolumns
csc                <- rep(hmcol[50], 30)
csc[group == 'T']  <- hmcol[200]
# column side color will be purple for T and orange for B
heatmap(sig.eset, scale = "row", col = hmcol, ColSideColors = csc)


#------------------------------------- w5 → w6


#LDA

#install.packages("BiocManager")
#BiocManager::install("genefilter")
library("genefilter")
#BiocManager::install("GEOquery")
library("GEOquery")

gse <- getGEO("GSE41526")

length(gse)
gse <- gse[[1]]
show(gse)

ex <- exprs(gse)
dim(ex)
colnames(ex)

# ex2 <- log2(ex)
#ex2 <- log2(na.omit(ex) + 40) #da non fare senza buon motivo
ex2 <- na.omit(ex)
ex3 <- ex2[ , 1:40] # prendiamo 40 campioni perhé sno quelli relativi al cancro al seno. Tutti gli altri sono per diverse malattie diversi

# sappiamo che sono così divisi nel dataset originale
f <- factor(c(rep(0, 20), rep(1, 20)))

tt40 <- rowttests(ex3, f) #dal package genefilter
#fa un test riga per riga

#tengo solo le righe che hanno p-value come indicato
keepers <- which(tt40$p.value < 0.1)
#keepers <- which(p.adjust(tt40$p.value) < 0.1)
#p.adjust() a cosa serve? se non lo si fa, si ottengono molte più features
#se lo si usa aumentare la soglia perché altrimenti non si hanno geni da utilizzare

tex3 <- t(ex3)
tex3 <- tex3[ , keepers]

dat <- cbind(tex3, c(rep(0, 20), rep(1, 20)))
colnames(dat)[ncol(dat)] <- "AFFECTED"

n.controls <- 20
n.affected <- 20
# IMPORTANTE garantire lo stesso numero
train <- sample( 1:(n.controls), (n.controls - 5))
test  <- setdiff(1:(n.controls), train)

test  <- c(test , test  + 20)
train <- c(train, train + 20)

library("MASS")
mod <- lda(AFFECTED ~ ., data = as.data.frame(dat), 
		   prior = c(0.5, 0.5), subset = train)
which(tt40$p.value < 0.1)
plot(mod)

mod.values <- predict(mod, as.data.frame(dat[train, ]))

plot(mod.values$x[ , 1], ylab = c("LDA Axis"))
text(mod.values$x[ , 1], col  = c(as.numeric(dat[train, "AFFECTED"]) + 10))

preds <- predict(mod, as.data.frame(dat[test, ]))
preds$class

table(as.numeric(preds$class),
	  as.numeric(dat[test,"AFFECTED"]))



#ROC CURVE
#install.packages("BiocManager")
#BiocManager::install("genefilter")
library("genefilter")
#BiocManager::install("GEOquery")
library("GEOquery")

library("MASS")
mod <- lda(AFFECTED ~ ., data = as.data.frame(dat),
		   prior = c(0.5, 0.5), subset = train)

plot(mod)

mod.values <- predict(mod, as.data.frame(dat[train, ]))

mod.values$class

plot(mod.values$x[ , 1], ylab=c("LDA Axis"))
text(mod.values$x[ , 1], col = c(as.numeric(dat[train, "AFFECTED"]) + 10))

preds <- predict(mod, as.data.frame(dat[test, ]))
preds$class

table(as.numeric(preds$class),
	  as.numeric(dat[test, AFFECTED]))





#BiocManager::install("pROC")
library("pROC")

roc_lda <- plot.roc(as.numeric(preds$class), 
					as.numeric(dat[test, "AFFECTED"]) )
plot(roc_lda, col = "grey")


#CARET
#install.packages("BiocManager")
#BiocManager::install("genefilter")
library("genefilter")
#BiocManager::install("GEOquery")
library("GEOquery")

gse <- getGEO("GSE41526")

length(gse)
gse <- gse[[1]]
show(gse)

ex <- exprs(gse)
dim(ex)
colnames(ex)

# ex2 <- log2(ex)
ex2 <- log2(na.omit(ex) + 40)

ex3 <- ex2[ , 1:40]

library("genefilter")
f       <- factor(c(rep(0, 20), rep(1, 20)))
tt40    <- rowttests(ex3, f)
keepers <- which(p.adjust(tt40$p.value) < 0.1)

tex3 <- t(ex3)

tex3 <- tex3[ , keepers]

dat  <- cbind(tex3, c(rep(0, 20), rep(1, 20)))
colnames(dat)[ncol(dat)] <- "AFFECTED"

#install.packages("caret")
library("caret")

inTrain <- createDataPartition(y = dat, p = .75, list = FALSE)

train <- dat[ inTrain, ]
test  <- dat[-inTrain, ]

Ctrl <- trainControl(method     = "repeatedcv",
					 numbers    = 4           , 
					 repeats    = 3           ,
					 classProbs = TRUE        ,
					 summaryFunction = twoClassSummary)

ldaFit <- train(AFFECTED ~ .,
				data      = train,
				method    = "lda",
				trControl = ctrl,
				metric    = "ROC")

plot(ldaFit)

pred <- predict(ldaFit, newdata = test)