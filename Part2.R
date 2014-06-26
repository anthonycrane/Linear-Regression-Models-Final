# Set seed/load packages
set.seed(42)
library(quantreg)
library(scales)
library(mgcv)


# Generate data with normal errors and introduce an outlier to observe 
# change in beta estimates in OLS and quantile regression              

X <- sort(runif(50,0,10))
eps <- rnorm(50)
Y <- 2*X + eps

linmod <- lm(Y~X)
qmod <- rq(Y~X)

plot(Y~X, pch=19)
abline(linmod,lty=2, col='red')
abline(qmod,lty=1,col='red')

X <- c(X,15)
Y <- c(Y,3)

linmod.out <- lm(Y~X)
qmod.out <- rq(Y~X)

plot(Y~X, pch=19)
abline(linmod.out,lty=2, col='red')
abline(qmod.out,lty=1,col='red')

summary(linmod)
summary(linmod.out)

summary(qmod)
summary(qmod.out)


# Run 10000 regressions of median and OLS, save beta estimates and plot them

betas <- matrix(ncol=2,nrow=10000)

for(i in 1:10000){
  X <- sort(runif(100,0,10))
  eps <- rnorm(100)
  Y <- 2*X + eps
  
  linmod <- lm(Y~X)
  qmod <- rq(Y~X)
  
  betas[i,1] <- coef(linmod)[[2]]
  betas[i,2] <- coef(qmod)[[2]]
}

plot(density(betas[,1]),xlim=c(1.8,2.2),xlab=expression(widehat(beta)[1]), main='',col='red')
lines(density(betas[,2]))

# Run 10000 regression of median and OLS introducing a random outlier, save beta estimates and plot them

betas.out <- matrix(ncol=2,nrow=10000)

for(i in 1:10000){
  X <- sort(runif(100,0,10))
  eps <- rnorm(100)
  Y <- 2*X + eps
  
  X <- c(X,rnorm(1,15,3))
  Y <- c(Y,rnorm(1,3,1))
  
  linmod <- lm(Y~X)
  qmod <- rq(Y~X)
  
  betas.out[i,1] <- coef(linmod)[[2]]
  betas.out[i,2] <- coef(qmod)[[2]]
}

plot(density(betas.out[,2]),xlim=c(1.2,2.2),xlab=expression(widehat(beta)[1]), main='')
lines(density(betas.out[,1]),col='red')

# Fit gam with quantile regression
load(file='P8111_FinalProjectData.rda')

# Drop 'parity' because only 3 non-zero observations, other data cleaning explained in Part1.R
data <- data[-which(names(data) == 'parity')]
data <- data[-which(names(data) %in% c('ppwt','delwt'))]
data <- data[-which(data$menarche==0),]
data$babysex <- as.factor(data$babysex)
data$frace <- as.factor(data$frace)
data$malform <- as.factor(data$malform)
data$mrace <- as.factor(data$mrace)

# Fit the model in Part1 and fit the equivalent median regression model using rqss
model.add <- gam(bwt ~ babysex + bhead + blength + mheight + mrace + ppbmi + smoken + fincome + s(wtgain) + s(gaweeks), data=data)
model.qadd <- rqss(bwt ~ babysex + bhead + blength + mheight + mrace + ppbmi + smoken + fincome + qss(wtgain,lambda=8) + qss(gaweeks,lambda=8),data=data)

summary(model)

# Calculate and plot residuals of quantile additive models
res <- data$bwt - fitted(model.qadd)
scatter.smooth(res~fitted(model.qadd), pch = 19,lpars=c(col='red'),ylim=c(-1500,2500),xlim=c(-100,4500))
qqnorm(res, pch = 19)
qqline(res,col='red')

# Generate plot of coefficient values across multiple quantile regressions

params <- matrix(ncol=6,nrow=11)

model.gam <- gam(bwt ~ babysex + bhead + blength + mheight + mrace + ppbmi + smoken + fincome + s(wtgain) + s(gaweeks), data=data)
for(i in 1:11){
  params[i,1] <- model.gam$coef[[i]]
}

tau <- c(0.1,0.25,0.5,0.75,0.9)
k <- 1
for(i in tau){
  model.tau <- rqss(bwt ~ babysex + bhead + blength + mheight + mrace + ppbmi + smoken + fincome + qss(wtgain,lambda=8) + qss(gaweeks,lambda=8),tau=i,data=data)
  k <- k+1
  for(j in 1:11){
    params[j,k] <- model.tau$coef[[j]]
  }
}

CI.upper  <- matrix(ncol=5,nrow=11)
tau <- c(0.1,0.25,0.5,0.75,0.9)
k <- 0
for(i in tau){
  model.tau <- rqss(bwt ~ babysex + bhead + blength + mheight + mrace + ppbmi + smoken + fincome + qss(wtgain,lambda=8) + qss(gaweeks,lambda=8),tau=i,data=data)
  sum <- summary(model.tau)
  k <- k+1
  for(j in 1:11){
    CI.upper[j,k] <- sum$coef[j,1] + 1.96*sum$coef[j,2]
  }
}

CI.lower  <- matrix(ncol=5,nrow=11)
tau <- c(0.1,0.25,0.5,0.75,0.9)
k <- 0
for(i in tau){
  model.tau <- rqss(bwt ~ babysex + bhead + blength + mheight + mrace + ppbmi + smoken + fincome + qss(wtgain,lambda=8) + qss(gaweeks,lambda=8),tau=i,data=data)
  sum <- summary(model.tau)
  k <- k+1
  for(j in 1:11){
    CI.lower[j,k] <- sum$coef[j,1] - 1.96*sum$coef[j,2]
  }
}

rownames(params) <- c("(Intercept)", "Sex", "Head Length", "Baby Length","Mother Height", 
                      "Mother Race:\n Black", "Mother Race:\n Asian","Mother Race:\n Puerto Rican",
                      "Pre-Pregnancy BMI","Cigarettes/Day","Family Income")
colnames(params) <- c('Mean','10th','25th','50th','75th','90th')

model.gam <- gam(bwt ~ babysex + bhead + blength + mheight + mrace + ppbmi + smoken + fincome + s(wtgain) + s(gaweeks), data=data)

par(mfrow=c(2,5),mar=c(5,5,4,2),las=2)
for(i in 2:11){
  upper.gam <- summary(model.gam)$p.coeff[[i]] + 1.96*summary(model.gam)$se[[i]]
  lower.gam <- summary(model.gam)$p.coeff[[i]] - 1.96*summary(model.gam)$se[[i]]
  plot(params[i,2:6],type='n',main=paste(rownames(params)[i]),
       ylim=c(min(c(lower.gam,upper.gam,CI.lower[i,],CI.upper[i,])),max(c(lower.gam,upper.gam,CI.lower[i,],CI.upper[i,]))),
       xlab='Percentile',ylab=expression(widehat(beta)),axes=F)
  polygon(c(1:5, rev(1:5)), c(CI.upper[i,], rev(CI.lower[i,])),col = alpha('red',0.2), border = NA)
  points(params[i,2:6],pch=19,type='b')
  abline(h=model.gam$coefficients[i],lty=1)
  abline(h=upper.gam,lty=2)
  abline(h=lower.gam,lty=2)
  axis(1, at = 1:5, labels = colnames(params)[2:6], las = 1)
  axis(2)
  box()
}

dev.off()

# Compare models via cross validation
CV = matrix(nrow=100,ncol=2,data=rep(NA,500))
colnames(CV) <- c("Median","GAM")

i <- 1
while(i < 101){
  GROUP = sample(1:5, size = nrow(data), replace = TRUE)
  data.train = data[which(GROUP != 1),]
  data.test = data[which(GROUP == 1),]
  
  # Predict for rqss won't predict extrapolations; must loop over randomly generated datasets and try predict. If it fails, continue and try again
  tryCatch({
    med <- rqss(bwt ~ babysex + bhead + blength + mheight + mrace + ppbmi + smoken + fincome + qss(wtgain,lambda=8) + qss(gaweeks,lambda=8),data=data.train)
    CV[i,1] <- mean((data.test$bwt - predict(med, newdata = data.test[,-4]))^2)
  
    gam3 <- gam(bwt ~ babysex + bhead + blength + mheight + mrace + ppbmi + smoken + fincome + s(wtgain) + s(gaweeks), select=T,data=data.train)
    CV[i,2] <- mean((data.test$bwt - predict(gam3, newdata = data.test[,-4]))^2)
    
    i <- i + 1
  }, error=function(e){}
  )
}

plot(1:2, CV[1,], axes = FALSE, type = 'o', pch = 19, 
     ylim=c(min(CV),max(CV)), xlab = "",ylab="CV",col=alpha('black',0.2))
axis(1, at = 1:2, labels = colnames(CV), las = 2)
axis(2)
box()

for(i in 2:99){
  points(CV[i,],type='b',pch=19, col=alpha('black',0.2))
}
