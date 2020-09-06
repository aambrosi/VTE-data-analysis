#################################
# Data Analysis and Exploration #
# A.A. 2019-2020                #
# Ambrosi Andrea                #
# Venous Thromboembolism        #
#################################

library("GEOquery")
gse  <- getGEO(getGPL = FALSE, GEO = "GSE48000",  filename = "/home/andrea/Desktop/data_analysis/GSE48000_series_matrix.txt.gz")
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
head(ex, 5)
colnames(ex)
summary(ex)
#    Analyze value distributions
boxplot(ex) #e un casino da vedere, troppe colonne e troppi dati


## --------------------------------------------------------------------------


# PCA
pca <- prcomp(t(ex)) #observations as rows and features as columns
pca$sdev
summary(pca)
screeplot(pca)

var_explained_df <- data.frame(PC=(1:132),
							   var_explained=(pca$sdev)^2/sum((pca$sdev)^2))
var_explained_df$PC <- as.factor(var_explained_df$PC)
head(var_explained_df)
var_explained_df <- var_explained_df[1:10,]
library(dplyr)
var_explained_df %>%
	ggplot(aes(x = PC, y=var_explained))+
	geom_col()+
	labs(y="Variance explained (%)") + theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 0))
pca$x
# draw PCA plot
# grey = Healthy Control, red = High Risk, green = Low Risk, yellow = Moderate Risk
grpcol <- c(rep("grey", 25), rep("red", 40), rep("green", 34), rep("yellow", 33))
plot(pca$x[ , 1], pca$x[ , 2], xlab = "Principal component 1", ylab = "Principal component 2", 
	 main = "PCA 1 vs. 2", type = "p", pch = 10, col = grpcol)

library(ggplot2)

pcadata <- data.frame(pca$x, Group=gse$`phenotype/group:ch1`)
head(pcadata)
ggplot(pcadata, aes(x = PC1, y = PC2, fill = Group, col = Group)) + 
	geom_point() + 
	labs(x = "Principal component 1", 
		 y = "Principal component 2") + 
	scale_x_continuous(position = "top") + 
	theme(legend.position = "bottom", legend.title = element_blank(), ) + 
	scale_colour_manual(values = c("#F8766D", "#7CAE00", "#C77CFF", "#00BFC4"))
#text(pca$x[ , 1], pca$x[ , 2], rownames(pca$x), cex = 0.75)


## --------------------------------------------------------------------------


#R CODE –K-MEANS EXAMPLE
#BiocManager::install("useful")
library("useful")

k <- 4
set.seed(1)
kmeans_result <- kmeans(t(ex), k)

table(kmeans_result$cluster)

#cluster visualization reduction using PCA 

plot(kmeans_result, data = t(ex))

cl1 <- names(kmeans_result$cluster[kmeans_result$cluster == 1])
cl2 <- names(kmeans_result$cluster[kmeans_result$cluster == 2])
cl3 <- names(kmeans_result$cluster[kmeans_result$cluster == 3])
cl4 <- names(kmeans_result$cluster[kmeans_result$cluster == 4])

prova <- data.frame(rbind(cbind(1, cl1), cbind(2, cl2), cbind(3, cl3), cbind(4, cl4)))
prova$cl1 <- as.character(prova$cl1)
prova <- prova[order(prova$cl1), ]

pcadata <- pcadata[,c("PC1","PC2","Group")]
pcadata <- cbind(pcadata, "name" = rownames(pcadata))

prova <- merge(pcadata, prova, by.x = "name", by.y = "cl1", all.x = TRUE, all.y = TRUE)
colnames(prova)[5] <- "K means classification"
ggplot(prova) + 
	geom_point(aes(x = PC1, y = PC2, fill = `K means classification`, col = `K means classification`)) + 
	labs(x = "Principal component 1", 
		 y = "Principal component 2") + 
	theme(legend.position = "bottom", legend.title = element_blank()) + 
	scale_x_continuous(position = "top") + 
	scale_colour_manual(values = c("#F8766D", "#7CAE00", "#C77CFF", "#00BFC4"))

prova <- table(prova$`K means classification`, prova$Group) #funziona male ma comunque sembra che distingua bene in due casi
1 - sum(diag(prova))/sum(prova) # error rate

#R CODE –HIERARCHICAL CLUSTERING EXAMPLE

dist_matrix <- dist(t(ex))
hc_result_ward.D   <- hclust(dist_matrix, method = "ward.D")
hc_result_ward.D2  <- hclust(dist_matrix, method = "ward.D2")
hc_result_single   <- hclust(dist_matrix, method = "single")
hc_result_complete <- hclust(dist_matrix, method = "complete")
hc_result_average  <- hclust(dist_matrix, method = "average")
hc_result_mcquitty <- hclust(dist_matrix, method = "mcquitty")
hc_result_median   <- hclust(dist_matrix, method = "median")
hc_result_centroid <- hclust(dist_matrix, method = "centroid")
#provare metodi diversi e commentare
# the agglomeration method to be used. This should be 
#(an unambiguous abbreviation of) one of "ward.D", 
#"ward.D2", "single", "complete", "average" (= UPGMA), 
#"mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).

k <- 4
groups_ward.D   <- cutree(hc_result_ward.D, k = k)
groups_ward.D2  <- cutree(hc_result_ward.D2, k = k)
groups_single   <- cutree(hc_result_single, k = k)
groups_complete <- cutree(hc_result_complete, k = k)
groups_average  <- cutree(hc_result_average, k = k)
groups_mcquitty <- cutree(hc_result_mcquitty, k = k)
groups_median   <- cutree(hc_result_median, k = k)
groups_centroid <- cutree(hc_result_centroid, k = k)

hc_results <- data.frame(table(groups_ward.D))
colnames(hc_results)[2] <- "Ward D"
hc_results <- merge(hc_results, data.frame(table(groups_ward.D2)), 
					by.x = "groups_ward.D", by.y = "groups_ward.D2")
colnames(hc_results)[3] <- "Ward D2"
hc_results <- merge(hc_results, data.frame(table(groups_single)), 
					by.x = "groups_ward.D", by.y = "groups_single")
colnames(hc_results)[4] <- "Single"
hc_results <- merge(hc_results, data.frame(table(groups_complete)), 
					by.x = "groups_ward.D", by.y = "groups_complete")
colnames(hc_results)[5] <- "Complete"
hc_results <- merge(hc_results, data.frame(table(groups_average)), 
					by.x = "groups_ward.D", by.y = "groups_average")
colnames(hc_results)[6] <- "Average"
hc_results <- merge(hc_results, data.frame(table(groups_mcquitty)), 
					by.x = "groups_ward.D", by.y = "groups_mcquitty")
colnames(hc_results)[7] <- "Mcquitty"
hc_results <- merge(hc_results, data.frame(table(groups_median)), 
					by.x = "groups_ward.D", by.y = "groups_median")
colnames(hc_results)[8] <- "Median"
hc_results <- merge(hc_results, data.frame(table(groups_centroid)), 
					by.x = "groups_ward.D", by.y = "groups_centroid")
colnames(hc_results)[9] <- "Centroid"
colnames(hc_results)[1] <- "Group"
data.frame(t(hc_results))

plot(hc_result_ward.D, hang = -1, labels = groups_ward.D)
rect.hclust(hc_result_ward.D, k = 4, which = NULL, x = NULL, h = NULL, border = 2, cluster = NULL) # red boxes to show groups

#ROC CURVE
#BiocManager::install("pROC")
library("pROC")

groups_ward.D2[groups_ward.D2==1]<- "Healthy Control"
groups_ward.D2[groups_ward.D2==2]<- "Low Risk"
groups_ward.D2[groups_ward.D2==3]<- "Moderate Risk"
groups_ward.D2[groups_ward.D2==4]<- "High Risk"
multiclass.roc(as.numeric(factor(groups_ward.D2)), as.numeric(pcadata$Group))
1-sum(diag(table(groups_ward.D2,pcadata$Group)))/sum(table(groups_ward.D2,pcadata$Group))
groups_complete[groups_complete==1]<- "Healthy Control"
groups_complete[groups_complete==3]<- "Low Risk"
groups_complete[groups_complete==2]<- "Moderate Risk"
groups_complete[groups_complete==4]<- "High Risk"
multiclass.roc(as.numeric(factor(groups_complete)), as.numeric(pcadata$Group))
1-sum(diag(table(groups_complete,pcadata$Group)))/sum(table(groups_complete,pcadata$Group))
groups_mcquitty[groups_mcquitty==1]<- "Healthy Control"
groups_mcquitty[groups_mcquitty==2]<- "Low Risk"
groups_mcquitty[groups_mcquitty==3]<- "Moderate Risk"
groups_mcquitty[groups_mcquitty==2]<- "High Risk"
multiclass.roc(as.numeric(factor(groups_mcquitty)), as.numeric(pcadata$Group))
1-sum(diag(table(groups_mcquitty,pcadata$Group)))/sum(table(groups_mcquitty,pcadata$Group))

#plottare sulla pca il risultato di ward complete e mcquitty
prova <- data.frame(pcadata, factor(groups_ward.D2), factor(groups_complete), factor(groups_mcquitty))
colnames(prova)[5:7] <- c("Ward D2 clusters", "Complete clusters", "Mcquitty clusters")
head(prova)
str(prova)
ggplot(prova, aes(x = PC1, y = PC2, colour = `Mcquitty clusters`)) + 
	geom_point() + 
	labs(x = "Principal component 1", 
		 y = "Principal component 2") + 
	theme(legend.position = "bottom", legend.title = element_blank()) + 
	scale_x_continuous(position = "top") + 
	scale_colour_manual(values = c("#F8766D", "#00BFC4", "#7CAE00", "#C77CFF"))


## --------------------------------------------------------------------------


library(randomForest) 
rfdata <- cbind(data.frame(t(ex)), "Group" = as.factor(gse$`phenotype/group:ch1`))
print(date())
set.seed(3) 
rf <- randomForest(x = rfdata[ ,1:ncol(rfdata)-1], y = rfdata$Group, ntree = 1000)#, importance = TRUE)
print(date()) # it takes about 3 minutes
rf#too much errors ~56% with 1000 trees
#table(rf$predicted, rfdata$Group) 
importance       <- data.frame(rf$importance)
importance$names <- rownames(importance)
importance       <- importance[order(-importance$MeanDecreaseGini), ]

plot(importance$MeanDecreaseGini[1:1000])

write.csv(importance[1:400,]$names, "names2.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)

# the code below have been runned several times tuning the parameters.
# Now topPredictors <- importance[1:readings[i], ] is fixed to 400 given that is the
# best result achieved.

readings <- seq(1, 1000, by = 1)
res <- rep(0, length(readings))
roc <- rep(0, length(readings))
cl <- parallel::makePSOCKcluster(6)
doParallel::registerDoParallel(cl)
for(i in 1:length(readings)){
	set.seed(i)
	inTrain   <- caret::createDataPartition(groups, list = FALSE, p = 0.75)
	topPredictors <- importance[1:400, ]$names
	set.seed(i)
	rf <- randomForest(x     = rfdata[ inTrain, topPredictors], 
					   xtest = rfdata[-inTrain, topPredictors], 
					   y     = rfdata$Group[ inTrain], 
					   ytest = rfdata$Group[-inTrain],
					   ntree = 1000)
	res[i] <- 1 - (sum(diag(table(rf$test$predicted, rfdata$Group[-inTrain]))) / sum(table(rf$test$predicted, rfdata$Group[-inTrain])))
	roc[i] <- as.numeric(multiclass.roc(as.numeric(rf$test$predicted), as.numeric(rfdata$Group[-inTrain]))$auc)
	}
parallel::stopCluster(cl)
summary(res)
summary(roc)
rf$importance
#toPlot <- data.frame("predictors" = as.factor(readings), "seed1" = res)
#toPlot$seed1 <- res
#toPlot$seed3 <- res

#toPlot[toPlot$seed3==min(min(toPlot$seed1), min(toPlot$seed2), min(toPlot$seed3)),] #top

ggplot(toPlot) + 
	geom_point(aes(x = predictors, y = seed1, group = 1, col = "1")) + 
	geom_line( aes(x = predictors, y = seed1, group = 1, col = "1")) + 
	geom_point(aes(x = predictors, y = seed2, group = 2, col = "2")) + 
	geom_line( aes(x = predictors, y = seed2, group = 2, col = "2")) + 
	geom_point(aes(x = predictors, y = seed3, group = 3, col = "3")) + 
	geom_line( aes(x = predictors, y = seed3, group = 3, col = "3")) + 
	labs(x = "# of predictors", 
		 y = "Error rate") + 
	scale_x_discrete(position = "top") + 
	theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0), 
		  legend.position = "bottom", legend.title = element_blank())

rfdata$PC1 <- pcadata$PC1
rfdata$PC2 <- pcadata$PC2
rfdata$predicted <- rf$predicted
ggplot(rfdata, aes(x = PC1, y = PC2, colour = predicted)) + 
	geom_point() + 
	labs(x = "Principal component 1", 
		 y = "Principal component 2") + 
	theme(legend.position = "bottom", legend.title = element_blank()) + 
	scale_x_continuous(position = "top")


## --------------------------------------------------------------------------


#LDA
library("MASS")
n.healthy  <- 25
n.high     <- 40
n.low      <- 34
n.moderate <- 33
p          <- 0.75 #probability of taking part to the training set
# IMPORTANTE garantire lo stesso numero
errors <- rep(0, 100)
auc    <- rep(0, 100)
i = 1
for(i in 2:100){
	set.seed(i)
	train <- sample(1:n.healthy, (round(p*n.healthy)))
	train <- append(train, sample(1:n.high    , (round(p*n.high    ))) + n.healthy)
	train <- append(train, sample(1:n.low     , (round(p*n.low     ))) + n.healthy + n.high)
	train <- append(train, sample(1:n.moderate, (round(p*n.moderate))) + n.healthy + n.high + n.low)
	test  <- setdiff(1:(n.healthy + n.high + n.low + n.moderate), train)
	train <- sort(train)
	
	dat <- data.frame(cbind(t(ex)))
	dat <- dat[, importance[1:400, ]$names]
	dat <- data.frame(apply(dat, 2, as.numeric))
	dat <- cbind(dat, "Group" = as.character(pcadata$Group))
	
	#table(dat$Group)/sum(table(dat$Group)) # so see the distribution of groups (maybe it would be usful for the prior)
	
	mod <- lda(pcadata$Group ~ ., data = dat, subset = train)
	#plot(mod)
	mod.values <- predict(mod, as.data.frame(dat[train, ]))
	
	#plot(mod.values$x[ , 1], ylab = c("LDA Axis"), col  = c(as.numeric(dat[train, "Group"]) + 10))
	#text(mod.values$x[ , 1], col  = c(as.numeric(dat[train, "Group"]) + 10))
	
	preds <- predict(mod, as.data.frame(dat[test, ]))
	#table(preds$class, dat[test, "Group"])
	
	r <- table(preds$class, dat[test, "Group"])
	
	errors[i] <- 1 - sum(diag(r))/sum(r) # error rate
	auc[i]    <- as.numeric(multiclass.roc(as.numeric(preds$class), 
			  							   as.numeric(dat[test, "Group"]), col = i, add = TRUE)$auc)
}

#Best results with seed=97
set.seed(97)
train <- sample(1:n.healthy, (round(p*n.healthy)))
train <- append(train, sample(1:n.high    , (round(p*n.high    ))) + n.healthy)
train <- append(train, sample(1:n.low     , (round(p*n.low     ))) + n.healthy + n.high)
train <- append(train, sample(1:n.moderate, (round(p*n.moderate))) + n.healthy + n.high + n.low)
test  <- setdiff(1:(n.healthy + n.high + n.low + n.moderate), train)
train <- sort(train)
dat <- data.frame(cbind(t(ex)))
dat <- dat[, importance[1:400, ]$names]
dat <- data.frame(apply(dat, 2, as.numeric))
dat <- cbind(dat, "Group" = as.character(pcadata$Group))
mod <- lda(pcadata$Group ~ ., data = dat, subset = train)
mod.values <- predict(mod, as.data.frame(dat[train, ]))
preds <- predict(mod, as.data.frame(dat[test, ]))
table(preds$class, dat[test, "Group"])
r <- table(preds$class, dat[test, "Group"])
1 - sum(diag(r))/sum(r) # error rate
as.numeric(multiclass.roc(as.numeric(preds$class), as.numeric(dat[test, "Group"]), col = i, add = TRUE)$auc)

#CARET
library("caret")
library("e1071")
## Run algorithms using 10-fold cross validation#
control <- trainControl(method = "cv", number = 35)
metric  <- "Accuracy"

fit.lda <- train(Group ~ . , data = dat, method = "lda", metric = metric, trControl = control)
fit.rf  <- train(Group ~ . , data = dat, method = "rf" , metric = metric, trControl = control)
results <- resamples(list(LDA = fit.lda, RF = fit.rf))
summary(results)
ggplot( results) + labs(y = "Accuracy")  

##Run algorithms using 10-fold cross validation, 10 times#
control   <- trainControl(method = "repeatedcv", number = 10, repeats = 10)
fit.lda.2 <- train(Group ~ . , data = dat, method = "lda", metric = metric, trControl = control)
fit.rf.2  <- train(Group ~ . , data = dat, method = "rf" , metric = metric, trControl = control)
results   <- resamples(list(LDA = fit.lda.2, RF = fit.rf.2))
ggplot(results) + labs(y = "Accuracy")


## --------------------------------------------------------------------------


library(rScudo)
library(caret)

groups <- pcadata$Group
ex     <- ex
set.seed(1)
inTrain   <- caret::createDataPartition(groups, list = FALSE, p = 0.75)
trainData <- ex[,  inTrain]
testData  <- ex[, -inTrain]

virtControl   <- rowMeans(trainData)
trainDataNorm <- trainData / virtControl
pVals         <- apply(trainDataNorm, 1, 
					   function(x) {
					   		stats::kruskal.test(x, groups[inTrain])$p.value
					   	})
trainDataNorm <- t(trainDataNorm[pVals <= 0.01, ])

cl <- parallel::makePSOCKcluster(6)
doParallel::registerDoParallel(cl)

try     <- c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500)
model   <- scudoModel(nTop = try, nBottom = try, N = 0.25)
set.seed(1)
control <- caret::trainControl(method = "cv", number = 10, summaryFunction = caret::multiClassSummary)
cvRes   <- caret::train(x = trainDataNorm, y = groups[inTrain], method = model, trControl = control)

parallel::stopCluster(cl)

model <- scudoClassify(ex[, inTrain], ex[, -inTrain], 0.25,
					   cvRes$bestTune$nTop, cvRes$bestTune$nBottom, groups[inTrain], alpha = 0.01)
caret::confusionMatrix(model$predicted, groups[-inTrain])

table(model$predicted, groups[-inTrain])

trainRes <- scudoTrain(trainData, groups = groups[inTrain], nTop = 300,
					   nBottom = 500, alpha = 0.01)
trainRes

trainNet <- scudoNetwork(trainRes, N = 0.25)
set.seed(2)
scudoPlot(trainNet, vertex.label = NA, x = c(-1,-1), y = c(-1.2,-1.2))
set.seed(2)
trainClust <- igraph::cluster_spinglass(trainNet, spins = 4)
set.seed(2)
plot(trainClust, trainNet, vertex.label = NA)

testRes <- scudoTest(trainRes, testData, groups[-inTrain], nTop = cvRes$bestTune$nTop,
					 nBottom = cvRes$bestTune$nBottom)
testRes
testNet <- scudoNetwork(testRes, N = 0.25)
set.seed(2)
scudoPlot(testNet, vertex.label = NA, x = c(-1,-1), y = c(-1.2,-1.2))
set.seed(2)
testClust <- igraph::cluster_spinglass(testNet, spins = 4)
set.seed(2)
plot(testClust, testNet, vertex.label = NA)

upSignatures(trainRes)[1:5, 1:5]
write.csv(consensusUpSignatures(testRes), col.names = TRUE, quote = FALSE, file = "up.csv")
write.csv(consensusDownSignatures(testRes), col.names = TRUE, quote = FALSE, file = "down.csv")

rScudo::scudoParams(trainRes)

classRes <- scudoClassify(trainData, testData, N = 0.25, nTop = cvRes$bestTune$nTop,
						  nBottom = cvRes$bestTune$nBottom, trainGroups = groups[inTrain], alpha = 0.01)
caret::confusionMatrix(classRes$predicted, groups[-inTrain])
multiclass.roc(as.numeric(model$predicted), as.numeric(groups[-inTrain]))
