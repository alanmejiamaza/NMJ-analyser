


DATABASE_full <-read.csv("/Users/alanmejiamaza/Desktop/Volocity NMJ analyser database.csv", header = TRUE)

DATABASE <-read.csv("/Users/alanmejiamaza/Desktop/DATABASE FINAL COMPLETE.csv", header = TRUE)

View(DATABASE)

install.packages("prcomp")
install.packages(c("FactoMineR", "factoextra"))
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/factoextra")
install.packages("caret", dependencies=c("Depends", "Suggests"))
install.packages("ggsignif")
DATABASE NAME: average_Gars_mice
GARS_ML
install.packages("imputePCA")
install.packages("prcomp")

1. # INTRODUCING THE DATABASE
library(ggplot2)
library(factoextra)
library(FactoMineR)
library("imputePCA")
library("prcomp")
library(PCAtools)
library(ggrepel)
library(reshape2)
library(lattice)
library(ggplot2)
library(factoextra)
library(FactoMineR)
library(BiocGenerics)
library(parallel)
library(corrplot)
library("ggpubr")
library("ggplot2")
library("magrittr")
library("ggpubr")
attach("BiocGenerics")
## nb <- estim_ncpPCA(orange,ncp.max=5) ## Time consuming, nb = 2

prcomp(na.omit(Volocity.database), center = TRUE, scale = TRUE)
prcomp(~V1+V2, data=Volocity.database, center = TRUE, scale = TRUE, na.action = na.omit)



View(Volocity.database)
Gars_active <- c(ALLDATABASE_PCA[,4:14])
View(Gars_active)

sum(is.na(PCA1.3m$genotype))
sum(is.na(PCA1.3m$timepoint))
sum(is.na(PCA1.3m$NMJ.counting))
sum(is.na(PCA1.3m$volume.red))
sum(is.na(PCA1.3m$volume.green))
sum(is.na(PCA1.3m$coverage))
sum(is.na(PCA1.3m$surface.red))
sum(is.na(PCA1.3m$surface.green))
sum(is.na(PCA1.3m$shape.factor.red))
sum(is.na(PCA1.3m$shape.factor.green))
sum(is.na(PCA1.3m$skeletal.length.red))
sum(is.na(PCA1.3m$skeletal.length.green))
sum(is.na(PCA1.3m$skeletal.diameter.red))
sum(is.na(PCA1.3m$skeletal.diameter.green))
na_count <- data.frame(na_count)

View(PCA1.3m)
View(PCA_FINAL_ALLDATABASE)
COMPLETE.DATABASE
COMPLETE.DATABASE <- na.omit(DATABASE) 

Gars_active <- c(PCA1.3m[,4:14])
Gars_active <- c(COMPLETE.DATABASE[,-c(1:3)])
Gars_active <- c(COMPLETE.DATABASE[,-c(1:4)])

COMPLETE.DATABASE <- na.omit(DATABASE[1:5053,])     # all 1-3.5 moths 
COMPLETE.DATABASE <- na.omit(DATABASE[5054:5856,])  # fus and tdp43 , 12months

COMPLETE.DATABASE <- na.omit(DATABASE[c(1:1841,2482:3462,5279:5790,4839:5053,5053:5184, 4001:4408),]) 
#just for mutants tdp43:5054-5184, gars=4001:4408, 
# fus-3m=4839:5053, fus-12m=5279:5790, sod-1m=1:1841, sod-3.5m=2427:3462
COMPLETE.DATABASE <- na.omit(DATABASE[2482:3462,]) # sod3.5 m
COMPLETE.DATABASE <- na.omit(DATABASE[5054:5184,]) # tdp43:5054-5184
COMPLETE.DATABASE <- na.omit(DATABASE[5279:5790,]) # fus-12m=5279:5790
COMPLETE.DATABASE <- na.omit(DATABASE[4001:4408,]) # fGARS 4001:4408
COMPLETE.DATABASE <- na.omit(DATABASE[4001:4408,]) # FUS-3M = 5415:5672
# sod-3.5m=2427:3462   FUS-3M = 5415:5672
COMPLETE.DATABASE <- na.omit(DATABASE[2427:3462,]) #  sod-3.5m=2427:3462 
COMPLETE.DATABASE <- na.omit(DATABASE[5415:5672,]) # FUS-3M = 5415:5672
COMPLETE.DATABASE <- na.omit(DATABASE[c(2427:3462,5415:5671),]) # FUS-3M = 5415:5671 + sod1-3.5 month
Gars_active
COMPLETE.DATABASE
View(DATABASE)
res.pca <- PCA(Gars_active, graph = FALSE)
print(res.pca)
eig.val <- get_eigenvalue(res.pca)
eig.val
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 100))
var <- get_pca(res.pca)
var
fviz_pca_var(res.pca, col.var ="black")
# corrplot
corrplot(var$cos2, is.corr = FALSE) 
#cos2 plot
fviz_cos2(res.pca, choice = "var", axes = 1:2)
#plot variables pCA
library("corrplot")
corrplot(var$contrib, is.corr=FALSE) 


fviz_pca_var(res.pca,
             col.var = "cos2", #color by contribution to the PC
             gradient.cols = c("#00AFBB", "#E7B600", "#FC4E07"),
             repel = TRUE   # avoid overlapping
) + theme_classic()
fviz_pca_var(res.pca, alpha.var = "cos2")
corrplot(var$contrib, is.corr = FALSE)
fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B600", "#FC4E07")
)+ theme_classic()

#plots of quality and contribution
fviz_pca_ind(res.pca)
fviz_pca_ind(res.pca, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B600", "#FC4E07"),
             repel = TRUE   # avoid text overlapping
) + theme_classic()

iris.pca <- iris[,-5]
fviz_pca_ind(iris.pca, 
             geom.ind = "point", #show points only (no text)
             col.ind = iris$Species, #color by group
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE,   #concentration ellispes
             legend.title = "Groups"
) + theme_classic()





DATABASE.pca <- PCA(DATABASE[c(2427:3462,5415:5671), -c(1:4)], graph = FALSE) # FUS-3M = 5415:5672
Gars_active

fviz_pca_ind(DATABASE.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = DATABASE[c(2427:3462,5415:5671),]$NMJ.counting, # color by groups
             palette = c("#FC4E07", "#E7B800", "#0033FF", "#0033FF","#00AFBB"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "NMJ status",
             axis.title=element_text(size=16, face="bold"),
             ggtheme=theme(axis.text=element_text(size=14)),
)








View(DATABASE[,25:57])
iris.pca <- PCA(iris[,-5], graph = FALSE)
fviz_pca_ind(iris.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = iris$Species, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
newdata <- na.omit(mydata)
NEWDATA <- na.omit(PCA_FINAL_ALLDATABASE)

PCA_FINAL_ALLDATABASE.pca <- PCA(PCA_FINAL_ALLDATABASE[, 4:14], graph = FALSE)
fviz_pca_ind(PCA_FINAL_ALLDATABASE.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = PCA_FINAL_ALLDATABASE[,]$timepoint, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "NMJ status",
             axis.title=element_text(size=16, face="bold"),
             ggtheme=theme(axis.text=element_text(size=14)),
             
)
View(PCA.12m)
PCA_12month.pca <- PCA(PCA.12m[1:625, 4:14], graph = FALSE)
fviz_pca_ind(PCA_12month.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = PCA.12m[1:625,]$NMJ.counting, # color by groups
             palette = c("#FC4E07", "#E7B800", "#0033FF", "#0033FF"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "NMJ status",
             axis.title=element_text(size=16, face="bold"),
             ggtheme=theme(axis.text=element_text(size=14)),
             
)
View(PCA1.3m)
PCA1.3m.pca <- PCA(PCA1.3m[341:4864, 4:14], graph = FALSE)
fviz_pca_ind(PCA1.3m.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = PCA1.3m[341:4864,]$genotype, # color by groups
             palette = c("#FC4E07", "#E7B800", "#0033FF", "#0033FF", "#666666"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "NMJ status",
             axis.title=element_text(size=16, face="bold"),
             ggtheme=theme(axis.text=element_text(size=14)),
             
)

which (is.na(DATABASE$genotype))
which (is.na(DATABASE$timepoint))
which (is.na(DATABASE$NMJ.counting))
which (is.na(DATABASE$bounds_0red))
which (is.na(DATABASE$genotype))
which (is.na(DATABASE$genotype))
which (is.na(DATABASE$genotype))
which (is.na(DATABASE$genotype))


COMPLETE.DATABASE <- na.omit(DATABASE) 

COMPLETE.DATABASE <- na.omit(DATABASE[5054:5856,])  # fus and tdp43 , 12months
COMPLETE.DATABASE <- na.omit(DATABASE[1:5053,])     # all 1-3.5 moths 
COMPLETE.DATABASE <- na.omit(DATABASE[c(2427:3462, 4797:5053),]) #
COMPLETE.DATABASE <- na.omit(DATABASE[c(1:2424,4001:4796), ])  #1:2425, sod1, 1 month, gars 4001:4796

COMPLETE.DATABASE <- na.omit(DATABASE[c(1:1841,2482:3462,5279:5790,4839:5053,5053:5184, 4001:4408),]) 
                                #just for mutants tdp43:5054-5184, gars=4001:4408, 
                                 # fus-3m=4839:5053, fus-12m=5279:5790, sod-1m=1:1841, sod-3.5m=2482:3462


levels(DATABASE[1:5053,]$NMJ.counting)
percentage <- prop.table(table(DATABASE$NMJ.counting)) * 100 
cbind(freq=table(DATABASE$NMJ.counting), percentage=percentage)

" FOR DATABASE "
COMPLETE.DATABASE <- DATABASE
COMPLETE.DATABASE[,-c(1:3,40,56,59)]
View(COMPLETE.DATABASE)
View(DATABASE)
COMPLETE.DATABASE <- na.omit(DATABASE) 

COMPLETE.DATABASE.pca <- PCA(COMPLETE.DATABASE[,-c(1:3)], graph = FALSE)
fviz_pca_ind(COMPLETE.DATABASE.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = COMPLETE.DATABASE[,]$NMJ.counting, # color by groups
             palette = c("#FC4E07", "#E7B800", "#0033FF", "#003333", "#CCFF00"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "NMJ status",
             axis.title=element_text(size=16, face="bold"),
             ggtheme=theme(axis.text=element_text(size=14)),
) + theme_classic()












 2. #EIGENVALUES/VARIANCES
#FactoMineR scales the variables
PCA(x, scale.unit =TRUE, ncp=5, graph = TRUE)
View(Gars_active)
res.pca <- PCA(Volocity.database, graph = FALSE)
print(Volocity.database)
eig.val <-get_eigenvalue(Volocity.database)
eig.val
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50)) +  theme_classic()


################################
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('PCAtools')

install.packages("lattice")

devtools::install_github('kevinblighe/PCAtools')

PCA1.3m <-read.csv("/Users/alanmejiamaza/Desktop/final_plots/PCA/PCA1-3m.csv", header = TRUE)

which(is.na(Volocity.database$volume.red))
which(is.na(Volocity.database$volume.green))
which(is.na(Volocity.database$coverage))
which(is.na(Volocity.database$surface.red))
which(is.na(Volocity.database$surface.green))
which(is.na(Volocity.database$shape.factor.red))
which(is.na(Volocity.database$shape.factor.green))
which(is.na(Volocity.database$skeletal.length.red))
which(is.na(Volocity.database$skeletal.length.green))
which(is.na(Volocity.database$skeletal.diameter.red))
which(is.na(Volocity.database$skeletal.diameter.green))
which(is.na(Volocity.database$timepoint))
which(is.na(Volocity.database$NMJ.counting))
which(is.na(Volocity.database$genotype))

##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
# Machine learning 
View(PCA_FINAL_ALLDATABASE)
Data ("PCA_FINAL_ALLDATABASE")
colnames(SOD1_1M)
install.packages("caret")
install.packages("generics")
install.packages("gower")
library("caret")
dataset <- PCA_FINAL_ALLDATABASE [,1:14]
dataset <- PCA_FINAL_ALLDATABASE

View(dataset)
attach(dataset)
#create a list of 80% of the rows in the original dataset we can use for training
validation_index <- createDataPartition(dataset$NMJ.counting, p=0.8, list=FALSE)
# select 20% of the data for validation
validation <- dataset[-validation_index,]
# use the reamining 80% of data to training and testing the models
dataset <- dataset[validation_index,]
#dimensions of the dataset
dim(dataset)
#list types for each attribute 
sapply(dataset, class)
# take a peek at the first 5 rows of the data 
head(dataset)
#list of the levels for the class
levels(dataset$NMJ.counting)
#summarize the class distribution
percentage <- prop.table(table(dataset$NMJ.counting)) * 100 
cbind(freq=table(dataset$NMJ.counting), percentage=percentage)
#summarize attribute distributions
summary(dataset)
dev.off()
graphics.off()

#split input and output 
x <- dataset[,4:14]
y <- dataset[,3]
#boxplot for each attribute on one image
par(mfrow=c(4,14))
for (i in 4:14) {
  boxplot(x[,i], main=names(iris)[i])
}
#barplot for class breakdown 
plot(y)
# MULTIVARIABLE PLOTS
# scatterplot matrix
dev.off()

featurePlot(x=x, y=y, plot="ellipse")
#box and whisker plots for each attribute
featurePlot(x=x, y=y, plot= "box")
#density plots for each attribute bu class value
scales <- list(x=list(relation="free"), y=list(relation="free"))
featurePlot(x=x, y=y, plot = "density", scales=scales)


/////////////////////////////////////////
  # test harness
  #RUN ALGORITHMS USING 10-FOLD CROSS VALIDATION
control <- trainControl(method = "CV", number = 100)
metric <- "Accuracy"  
#TESTING 5 DIFFENRENT ALGORITHMS  
#LINEAR DISCRIMINANT ANALYSIS (LDA), LINEAR
#CLASSIFICATION AND REGRESSION TRESS(CART), NO LINEAR
#NEAREST NEIGHBORS(kNN), NO LINEAR
#SUPPORT VECTOR MACHINES(SVM) WITH A LINEAR KERNEL: COMPLEX NONLINEAR METHODS
#RANDOM FOREST (RF): COMPLEX NONLINEAR METHODS
dat <- na.omit(dat)
dataset <- na.omit(dataset)

apply(dataset, 2, function(x) any(is.na(x)))
ls()


# a) lenar algorithms
set.seed(7)
fit.lda <- train(NMJ.counting~., data=dataset, method="lda", metric=metric, trControl=control)
# b.1) nonlinear algorithms, CART
set.seed(7)
fit.cart <- train(NMJ.counting~., data=dataset, method="rpart", metric=metric, trControl=control)
# b.2) nonlinear algorithms, jNN
set.seed(7)
fit.knn <- train(NMJ.counting~., data=dataset, method="knn", metric=metric, trControl=control)
# c) advanced algorithms, SVM
set.seed(7)
fit.svm <- train(NMJ.counting~., data=dataset, method="svmRadial", metric=metric, trControl=control)
# d) random forest
set.seed(7)
fit.rf <- train(NMJ.counting~., data=dataset, method="rf", metric=metric, trControl=control)

# SUMMARIZA ACCURACY OF MODELS 
results <- resamples(list(lda=fit.lda, cart=fit.cart, knn=fit.knn, svm=fit.svm, rf=fit.rf))
summary(results)  
# compare accuracy of models
dotplot(results)
# summarize best model
print(fit.svm)

#estimate skill of LDA on the validation dataset
predictions <- predict(fit.rf, validation)
confusionMatrix(predictions, validation$NMJ.counting)

















