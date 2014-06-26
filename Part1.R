# Load required packages
library(glmnet)
library(mgcv)
library(scales)
library(ggplot2)

# Set seed and load data
set.seed(42)
setwd("~/Dropbox/regression/Final")
load(file='P8111_FinalProjectData.rda')

#########################
# Data Summary/Cleaning #
#########################
str(data)
summary(data)

# Drop variables 'pnumlbw' and 'pnumsga' as they do not have any information (all values are 0)
# Drop 'parity' because only 3 non-zero observations
# Drop 'delwt' and 'ppwt' as the difference between the two is coded in weight gain
# Remove subject who has age of menarche=0 (not biologically possible)

data <- data[-which(names(data) %in% c('pnumlbw','pnumsga'))]
data <- data[-which(names(data) == 'parity')]
data <- data[-which(names(data) %in% c('ppwt','delwt'))]
data <- data[-which(data$menarche==0),]

# Recode 'babysex', 'frace', 'malform', and 'mrace' as factors
data$babysex <- as.factor(data$babysex)
data$frace <- as.factor(data$frace)
data$malform <- as.factor(data$malform)
data$mrace <- as.factor(data$mrace)

########################
# Simple Linear Models #
########################

lm1 <- lm(data$bwt ~ ., data = data)
lm2 <- lm(data$bwt ~  . - ppbmi + poly(ppbmi,2), data = data)
lm3 <- lm(data$bwt ~  . - ppbmi + poly(ppbmi,3), data = data)
lm3 <- lm(data$bwt ~  . - wtgain + poly(wtgain,3), data = data)
lm4 <- lm(data$bwt ~  . - wtgain + poly(wtgain,3) - momage + poly(momage,2), data = data)
lm4 <- lm(data$bwt ~  . - wtgain + poly(wtgain,3) - momage + poly(momage,2) - menarche + poly(menarche,3), data = data)


################
# LASSO Models #
################

data.lasso <- data
data.lasso$frace2 <- ifelse(data.lasso$frace==2,1,0)
data.lasso$frace3 <- ifelse(data.lasso$frace==3,1,0)
data.lasso$frace4 <- ifelse(data.lasso$frace==4,1,0)
data.lasso$frace8 <- ifelse(data.lasso$frace==8,1,0)

data.lasso <- data.lasso[-which(names(data.lasso)=='frace')]

data.lasso$mrace2 <- ifelse(data.lasso$mrace==2,1,0)
data.lasso$mrace3 <- ifelse(data.lasso$mrace==3,1,0)
data.lasso$mrace4 <- ifelse(data.lasso$mrace==4,1,0)

data.lasso <- data.lasso[-which(names(data.lasso) %in% c('mrace','ppbmi'))]
data.lasso <- cbind(data.lasso,poly(data$ppbmi,3))
X <- data.matrix(data.lasso[,-4])
Y <- data.lasso[,4]

set.seed(1234)
CV <- cv.glmnet(X,Y)
CV$lambda.min

model.lasso1 <- glmnet(X,Y,lambda=1.44)
coef(model.lasso1)

##############
# GAM Models #
##############
model.1 <- gam(bwt ~ babysex + bhead + blength + frace  + mrace + mheight + ppbmi + smoken + fincome + s(wtgain) + s(gaweeks) + s(menarche), data=data)
AIC(model.1)
model.2 <- gam(bwt ~  malform + babysex + bhead + blength + mheight + mrace + ppbmi + smoken + fincome + s(wtgain) + s(gaweeks) + s(menarche), data=data)
AIC(model.2)
model.3 <- gam(bwt ~ babysex + bhead + blength + mheight + menarche + mrace + ppbmi + smoken + fincome + s(wtgain) + s(gaweeks), data=data)
AIC(model.3)
model.4 <- gam(bwt ~ babysex + bhead + blength + mheight + mrace + ppbmi + smoken + fincome + s(wtgain) + s(gaweeks) + s(momage), data=data)
AIC(model.4)
model.5 <- gam(bwt ~ babysex + bhead + blength + mheight + mrace + ppbmi + smoken + fincome + s(wtgain) + s(gaweeks), data=data)
AIC(model.5)



####################################
# Cross Validation of Various Fits #
####################################

CV = matrix(nrow=100,ncol=6,data=rep(NA,500))
colnames(CV) <- c("MLR1","MLR2", "LASSO", "GAM1","GAM2","GAM3")

for(i in 1:100){
  GROUP = sample(1:5, size = nrow(data), replace = TRUE)
  data.train = data[which(GROUP != 1),]
  data.test = data[which(GROUP == 1),]
  
  mlr1 <- lm(bwt ~ ., data=data.train)
  CV[i,1] <- mean((data.test$bwt - predict(mlr1, newdata = data.test[,-4]))^2)
  
  mlr2 <- lm(bwt ~ babysex + bhead + blength + mheight + mrace + ppbmi + smoken + fincome + poly(wtgain,3) + poly(gaweeks,3), data = data.train)  
  CV[i,2] <- mean((data.test$bwt - predict(mlr2, newdata = data.test[,-4]))^2)
  
  data.lasso <- data.train
  data.lasso$frace2 <- ifelse(data.lasso$frace==2,1,0)
  data.lasso$frace3 <- ifelse(data.lasso$frace==3,1,0)
  data.lasso$frace4 <- ifelse(data.lasso$frace==4,1,0)
  data.lasso$frace8 <- ifelse(data.lasso$frace==8,1,0)
  
  data.lasso <- data.lasso[-which(names(data.lasso)=='frace')]
  
  data.lasso$mrace2 <- ifelse(data.lasso$mrace==2,1,0)
  data.lasso$mrace3 <- ifelse(data.lasso$mrace==3,1,0)
  data.lasso$mrace4 <- ifelse(data.lasso$mrace==4,1,0)
  
  data.lasso <- cbind(data.lasso,poly(data.lasso$wtgain,3),poly(data.lasso$gaweeks,3))
  data.lasso <- data.lasso[-which(names(data.lasso) %in% c('wtgain','gaweeks'))]
  
  X <- data.matrix(data.lasso[,-4])
  Y <- data.lasso[,4]
  
  data.lasso.test <- data.test
  data.lasso.test$frace2 <- ifelse(data.lasso.test$frace==2,1,0)
  data.lasso.test$frace3 <- ifelse(data.lasso.test$frace==3,1,0)
  data.lasso.test$frace4 <- ifelse(data.lasso.test$frace==4,1,0)
  data.lasso.test$frace8 <- ifelse(data.lasso.test$frace==8,1,0)
  
  data.lasso.test <- data.lasso.test[-which(names(data.lasso.test)=='frace')]
  
  data.lasso.test$mrace2 <- ifelse(data.lasso.test$mrace==2,1,0)
  data.lasso.test$mrace3 <- ifelse(data.lasso.test$mrace==3,1,0)
  data.lasso.test$mrace4 <- ifelse(data.lasso.test$mrace==4,1,0)
  
  
  data.lasso.test <- cbind(data.lasso.test,poly(data.lasso.test$wtgain,3),poly(data.lasso.test$gaweeks,3))
  data.lasso.test <- data.lasso.test[-which(names(data.lasso.test) %in% c('wtgain','gaweeks'))]
  
  X.test <- data.matrix(data.lasso.test[,-4])
  Y.test <- data.lasso.test[,4]
  
  lasso1 <- glmnet(X,Y,lambda=1.44)
  CV[i,3] <- mean((Y.test - predict(lasso1, newx = X.test))^2)
  
  gam1 <- gam(bwt ~ babysex + bhead + blength + mheight + menarche + mrace + ppbmi + smoken + fincome + s(wtgain) + s(gaweeks), data=data.train)
  CV[i,4] <- mean((data.test$bwt - predict(gam1, newdata = data.test[,-4]))^2)
  
  gam2 <- gam(bwt ~ babysex + bhead + blength + mheight + mrace + ppbmi + smoken + fincome + s(wtgain) + s(gaweeks) + s(momage), data=data.train)
  CV[i,5] <- mean((data.test$bwt - predict(gam2, newdata = data.test[,-4]))^2)
  
  gam3 <- gam(bwt ~ babysex + bhead + blength + mheight + mrace + ppbmi + smoken + fincome + s(wtgain) + s(gaweeks), select=T,data=data.train)
  CV[i,6] <- mean((data.test$bwt - predict(gam3, newdata = data.test[,-4]))^2)

}

plot(1:6, CV[1,], axes = FALSE, type = 'o', pch = 19, 
     ylim=c(min(CV),max(CV)), xlab = "",ylab="CV",col=alpha('black',0.2))
axis(1, at = 1:6, labels = colnames(CV), las = 2)
axis(2)
box()

for(i in 2:100){
  points(CV[i,],type='b',pch=19, col=alpha('black',0.2))
}

###############
# Final Model #
###############
model <- gam(bwt ~ babysex + bhead + blength + mheight + mrace + ppbmi + smoken + fincome + s(wtgain) + s(gaweeks), data=data)
vis.gam(model, view=c('wtgain','gaweeks'),type='response', plot.type='contour',main='',contour.col='black')
plot(model)
summary(model)
coef(model)

scatter.smooth(model$fitted.values, model$residuals, pch = 19,lpars=c(col='red'))

qqnorm(model$residuals, pch = 19)
qqline(model$residuals,col='red')

plot(cooks.distance(model), pch=19)

#################
#################
## END OF FILE ##
#################
#################
