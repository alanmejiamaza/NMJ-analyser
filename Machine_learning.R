
# Machine learning project

Wed 29 January 2020

install.packages("caret")
# if need other packages use the following:
install.packages("caret", dependencies=c("Depends","Suggests"))
install.packages("ellipse")

library(lattice)
library(ggplot2)
library(caret)
library("ellipse")

# loading database 

SOD1_1M <-read.csv("/Users/alanmejiamaza/Desktop/final_plots/sod1_1m/SOD1_1M_DATABASE_CLEAN.csv", header = TRUE)
View(SOD1_1M)
#attach the iris dataset to the environment
data ("SOD1_1M")
colnames(SOD1_1M)
COMPLETE.DATABASE <- na.omit(DATABASE) 

dataset <- COMPLETE.DATABASE

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
x <- dataset[,4:56]
y <- dataset[,3]
#boxplot for each attribute on one image
par(mfrow=c(4,56))
for (i in 4:56) {
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
control <- trainControl(method = "CV", number = 10)
metric <- "Accuracy"  
#TESTING 5 DIFFENRENT ALGORITHMS  
#LINEAR DISCRIMINANT ANALYSIS (LDA), LINEAR
#CLASSIFICATION AND REGRESSION TRESS(CART), NO LINEAR
#NEAREST NEIGHBORS(kNN), NO LINEAR
#SUPPORT VECTOR MACHINES(SVM) WITH A LINEAR KERNEL: COMPLEX NONLINEAR METHODS
#RANDOM FOREST (RF): COMPLEX NONLINEAR METHODS
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
predictions <- predict(fit.svm, validation)

confusionMatrix(predictions, validation$NMJ.counting)
