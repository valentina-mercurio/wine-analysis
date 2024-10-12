##### PROGETTO MS-SL A.A. 2023/2024
# fonte dei dati : https://github.com/YBI-Foundation/Dataset/blob/main/Wine.csv
######## Caricamento data set ##########
rm(list=ls())
library(readr)
data <- read_csv("Wine.csv")
data<-as.data.frame(data)
str(data)
attach(data)
#View(data)
data<-data[,-2]

#View(data)
data<-data[,-2]
data$class_label<-factor(data$class_label, levels=c(1,2,3))
str(data)
######## Multicollinearità ######
names(data)
Cor<-cor(data[,3:ncol(data)])

panel.corrgram <- function(x, y, z, subscripts, at, level = 0.9, label = FALSE, ...) {
  #require("ellipse", quietly = TRUE)
  x <- as.numeric(x)[subscripts]
  y <- as.numeric(y)[subscripts]
  z <- as.numeric(z)[subscripts]
  zcol <- level.colors(z, at = at, ...)
  for (i in seq(along = z)) {
    ell <- ellipse(z[i], level = level, npoints = 50, scale = c(.2, .2),
                   centre = c(x[i], y[i]))
    panel.polygon(ell, col = zcol[i], border = zcol[i], ...)
  }
  if (label) {
    panel.text(x = x, y = y, lab = 100 * round(z, 2), cex = 0.8,
               col = ifelse(z < 0, "white", "black"))
  }
}
library(lattice); library(ellipse)
##### Stampare grafico
print(levelplot(Cor, at = do.breaks(c(-1.01, 1.01), 20), 
                xlab = NULL, ylab = NULL, colorkey = list(space = "top"), 
                scales = list(x = list(rot = 90)) , panel = panel.corrgram, 
                label = TRUE ,aspect = "fill",main ="Tabella di correlazione"))

x<-as.matrix(data[,3:ncol(data)])
xt<-t(x)
matr<-xt %*% x
Deter<-det(matr)
Deter

X<-as.matrix(cbind(rep(1,nrow(data)), data[,3:ncol(data)]))
autval<-eigen(t(X) %*% X)

MinAut<-min(autval$values)
min(autval$values)

condition_number<-sqrt(max(autval$values)/min(autval$values))
condition_number

Mult<-cbind(Deter, MinAut, condition_number)
colnames(Mult)<-c("Determinante", "   Min Autovalore", "   Condition Number")
Mult


######## Creazione modello con tecniche inferenziali ####

coeffvar<-flavanoids/total_phenols
mean(coeffvar)
QtàTotFlavFen<-round(((coeffvar*flavanoids)+total_phenols),3)
type_wine<-(factor(data$class_label, levels=c(1,2,3)))

type_wine<-(factor(data$class_label, levels=c(1,2,3)))


dataset <- with(data, {
  type_wine<-(factor(data$class_label, levels=c(1,2,3)))
  data.frame(type_wine,
             alcohol, malic_acid,ash,alcalinity_of_ash, 
             magnesium, QtàTotFlavFen, nonflavanoid_phenols,
             proanthocyanins,color_intensity,hue,od280,proline)
})

str(dataset)

#View(dataset)
mod<-lm(alcohol~.,data=dataset)
summary(mod)

names(dataset)
Cordataset<-cor(dataset[,3:ncol(dataset)])

library(GGally)
ggpairs(dataset[,3:ncol(dataset)])

##### Stampare grafico
print(levelplot(Cordataset, at = do.breaks(c(-1.01, 1.01), 20), 
                xlab = NULL, ylab = NULL, colorkey = list(space = "top"), 
                scales = list(x = list(rot = 90)) , panel = panel.corrgram, 
                label = TRUE ,aspect = "fill",main ="Tabella di correlazione"))


# creazione modello

mod<-lm(alcohol~.+poly(QtàTotFlavFen,2)-QtàTotFlavFen,data=dataset)
summary(mod)
mod1<-lm(log(alcohol)~.,data=dataset)
summary(mod1)
mod2<-lm(alcohol~.+QtàTotFlavFen*color_intensity,data=dataset)
summary(mod2)
mod3<-lm(alcohol~type_wine+malic_acid+color_intensity+proline+log(QtàTotFlavFen),data=dataset)
summary(mod3)


modello_finale<-lm((alcohol)~type_wine+malic_acid+ash+color_intensity*hue-color_intensity-hue,data=dataset)
summary(modello_finale)

#barolo, grignolino, barbera
library(MASS)
library(car)
mod_null <- lm(alcohol ~ 1, data=dataset)
stepAIC(mod_null, scope=list(lower=mod_null,upper=modello_completo), direction="forward", data=dataset)


modello_completo<-lm((alcohol)~.,data=dataset)
summary(modello_completo)
modelloSAR <- step(modello_completo, direction = "backward")
summary(modelloSAR)

modello_finale<-lm((alcohol)~type_wine+color_intensity,data=dataset)
summary(modello_finale)

######## TECNICHE DI REGOLARIZZAZIONE ##############

set.seed(100)

index=sample(1:nrow(data), 0.7*nrow(data))

trainData=data[index, ]
testData=data[-index, ]

dim(trainData);dim(testData)

mtrainData<-lm(alcohol~., data=trainData)
summary(mtrainData)

# CROSS VALIDATION
library(glmnet)
#names(data)
xx<-data[, -2]
x<-as.matrix (xx); dim (x)
y<-data[,2]

qq<-seq(10,-2,length=100)
griglia=10^qq 
# Scelta miglior metodo di Regolarizzazione per Min MSE predittivo
# 1. RR
ridge.mod.ALL=glmnet(x,y,alpha=0,lambda=griglia)
dim(coef(ridge.mod.ALL))

plot(ridge.mod.ALL, main="Ridge Regression\n", xvar="lambda",label=TRUE)

cv.outK10=cv.glmnet(x,y,lambda=griglia, alpha=0) 
plot(cv.outK10, main="k-fold CV RR\n")

bestLambda<-cv.outK10$lambda.min 
bestLambda
cv.outK10$lambda.1se

# calcolo modello col min lambda RR
ridge.mod.kCV=glmnet(x,y,alpha=0,lambda=bestLambda)
coef(ridge.mod.kCV)[,1]

## b) Ridge regression con il metodo leave one out cross validation:
n <- nrow(data)
cv.outLOOCV=cv.glmnet(x,y,lambda=griglia, nfolds=n, grouped=FALSE, alpha=0)

plot(cv.outLOOCV, main="LOOCV RR\n")

bestLambdaLOOCV<-cv.outLOOCV$lambda.min 
bestLambdaLOOCV
cv.outLOOCV$lambda.1se


# stima modello finale
ridge.mod.LOOCV=glmnet(x,y,alpha=0, lambda=bestLambdaLOOCV)
coef(ridge.mod.LOOCV)[,1]

# LASSO
LASSO.mod.ALL=glmnet(x,y,alpha=1, lambda=griglia)

plot(LASSO.mod.ALL, main="LASSO\n", xvar="lambda",label=TRUE)
dim(coef(LASSO.mod.ALL))

## a) Lasso con il metodo k-fold cross validation
cv.outK10.LASSO=cv.glmnet(x,y,lambda=griglia, alpha=1)

plot(cv.outK10.LASSO, main="LASSO: k-fold CV \n")

bestLambda.LASSO<-cv.outK10.LASSO$lambda.min
bestLambda.LASSO

#Stima del modello di regressione lasso con il valore minimo di lambda: 
LASSO.mod.kCV=glmnet(x,y,alpha=1,lambda=bestLambda.LASSO)
coef(LASSO.mod.kCV)[,1]

# ## b) Lasso con il metodo leave one out cross validation 
cv.outLOOCV.LASSO=cv.glmnet(x,y,lambda=griglia,nfolds=n, grouped=FALSE, alpha=1)

plot(cv.outLOOCV.LASSO, main="LASSO: LOOCV LASSO\n")

bestLambdaLOOCV.LASSO<-cv.outLOOCV.LASSO$lambda.min
bestLambdaLOOCV.LASSO

#È possibile, infine, effettuare un confronto tra le stime dei parametri ottenute con la tecnica lasso e ridge regression.
# Modello scelto con tecniche inferenziali
modello<-lm(alcohol~., data=data)
summary(modello)

modello$coefficients
ConStimeLASSORR<-cbind(coef(LASSO.mod.kCV)[,1], coef(ridge.mod.kCV)[,1])
colnames(ConStimeLASSORR)<-c("LASSO","RR")
ConStimeLASSORR

OLSCOEFF<-as.matrix(modello$coefficients)
LASSOCOEFF<-as.matrix(coef(LASSO.mod.kCV)[,1])

max_length <- max(length(OLSCOEFF), length(LASSOCOEFF))
length(OLSCOEFF) <- max_length                      
length(LASSOCOEFF) <- max_length

LASSO_OLS<-cbind(OLSCOEFF, LASSOCOEFF)
LASSO_OLS[is.na(LASSO_OLS)]<-0
colnames(LASSO_OLS)<-c("OLS","LASSO")

LASSO_OLS_RR<-cbind(LASSO_OLS, coef(ridge.mod.kCV)[,1])
colnames(LASSO_OLS_RR)<-c("OLS","LASSO","RR")
LASSO_OLS_RR

# EN
EN.modes.ALL <- glmnet(x, y, lambda=griglia, alpha=.1)
plot(EN.modes.ALL, main="ELASTIC NET; alpha=0.1\n",xvar="lambda",label=TRUE)

cv.outK10.EN01=cv.glmnet(x,y,lambda=griglia, alpha=0.1) 
plot(cv.outK10.EN01, main="Elastic Net alpha=0.1: k-fold CV\n")
bestLambda.EN01<-cv.outK10.EN01$lambda.min 
bestLambda.EN01
EN01.mod.kCV=glmnet(x,y,alpha=0.1,lambda=bestLambda.EN01)
coef(EN01.mod.kCV)[,1]

EN.modes.ALL <- glmnet(x, y, lambda=griglia, alpha=.3)
plot(EN.modes.ALL, main="ELASTIC NET; alpha=0.3\n",xvar="lambda",label=TRUE)

cv.outK10.EN03=cv.glmnet(x,y,lambda=griglia, alpha=0.3) 
plot(cv.outK10.EN03, main="Elastic Net alpha=0.3: k-fold CV \n")
bestLambda.EN03<-cv.outK10.EN03$lambda.min 
bestLambda.EN03
EN03.mod.kCV=glmnet(x,y,alpha=0.3,lambda=bestLambda.EN03)
coef(EN03.mod.kCV)[,1]

EN.modes.ALL <- glmnet(x, y, lambda=griglia, alpha=.5)
plot(EN.modes.ALL, main="ELASTIC NET; alpha=0.5\n",xvar="lambda",label=TRUE)

cv.outK10.EN05=cv.glmnet(x,y,lambda=griglia, alpha=0.5) 
plot(cv.outK10.EN05, main="Elastic Net alpha=0.5: k-fold CV \n")
bestLambda.EN05<-cv.outK10.EN05$lambda.min 
bestLambda.EN05
EN05.mod.kCV=glmnet(x,y,alpha=0.5,lambda=bestLambda.EN05)
coef(EN05.mod.kCV)[,1]

par(mfrow=c(1,2))
EN.modes.ALL <- glmnet(x, y, lambda=griglia, alpha=.7)
plot(EN.modes.ALL, main="ELASTIC NET; alpha=0.7\n",xvar="lambda",label=TRUE)

cv.outK10.EN07=cv.glmnet(x,y,lambda=griglia, alpha=0.7) 
plot(cv.outK10.EN07, main="Elastic Net alpha=0.7: k-fold CV \n")
bestLambda.EN07<-cv.outK10.EN07$lambda.min
bestLambda.EN07
EN07.mod.kCV=glmnet(x,y,alpha=0.7,lambda=bestLambda.EN07)
coef(EN07.mod.kCV)[,1]

EN.modes.ALL <- glmnet(x, y, lambda=griglia, alpha=.9)
plot(EN.modes.ALL, main="ELASTIC NET; alpha=0.9\n",xvar="lambda",label=TRUE)

cv.outK10.EN09=cv.glmnet(x,y,lambda=griglia, alpha=0.9) 
plot(cv.outK10.EN09, main="Elastic Net alpha=0.9: k-fold CV \n")
bestLambda.EN09<-cv.outK10.EN09$lambda.min 
bestLambda.EN09
EN09.mod.kCV=glmnet(x,y,alpha=0.9,lambda=bestLambda.EN09)
coef(EN09.mod.kCV)[,1]



# SCELTA DEL MODELLO
Tabcoeff<-cbind(coef(LASSO.mod.kCV)[,1], coef(ridge.mod.kCV)[,1], coef(EN01.mod.kCV)[,1],
                coef(EN03.mod.kCV)[,1], coef(EN05.mod.kCV)[,1], coef(EN07.mod.kCV)[,1],
                coef(EN09.mod.kCV)[,1] )
colnames(Tabcoeff)<-c("LASSO", "RR","EN 01","EN 03","EN 05","EN 07","EN 09")

#Esempio: per la stima lasso estraiamo gli mse

cv.outK10.LASSO$cvm
mse.minLASSO<-cv.outK10.LASSO$cvm[cv.outK10.LASSO$lambda == cv.outK10.LASSO$lambda.min]
mse.minRR<-cv.outK10$cvm[cv.outK10$lambda == cv.outK10$lambda.min]
mse.minEN01<-cv.outK10.EN01$cvm[cv.outK10.EN01$lambda == cv.outK10.EN01$lambda.min]
mse.minEN03<-cv.outK10.EN03$cvm[cv.outK10.EN03$lambda == cv.outK10.EN03$lambda.min]
mse.minEN05<-cv.outK10.EN05$cvm[cv.outK10.EN05$lambda == cv.outK10.EN05$lambda.min]
mse.minEN07<-cv.outK10.EN07$cvm[cv.outK10.EN07$lambda == cv.outK10.EN07$lambda.min]
mse.minEN09<-cv.outK10.EN09$cvm[cv.outK10.EN09$lambda == cv.outK10.EN09$lambda.min]

rbind(mse.minLASSO, mse.minRR, mse.minEN01, mse.minEN03, mse.minEN05, mse.minEN07, mse.minEN09)
min(mse.minLASSO, mse.minRR, mse.minEN01, mse.minEN03, mse.minEN05, mse.minEN07, mse.minEN09)

######## ETEROSCHEDASCITITA'##########
# ANALISI GRAFICA
# modelloSAR modello Selezione Automatica Regressori
res1<-resid(modelloSAR) ; fit1<-fitted(modelloSAR); plot(fit1,res1, main="Plot dei residui col modello scelto", ylab="Residui", xlab="Valori osservati")

#TEST BREUSCH-PAGAN
resS<-modelloSAR$residuals
resS2<-resS^2
modres12<-lm(resS2~., data=dataset)
summary(modres12)


# TEST DI WHITE
fit1<-fitted(modelloSAR)
fit12<-fit1**2
modresW<-lm(resS2~fit1+fit12)
summary(modresW)
