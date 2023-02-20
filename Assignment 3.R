library(mvtnorm)
library(mcmc)
library(tidyverse)

# Setting Seed
set.seed(1234)

#Reading the csv File
wine<- read.csv(("winequality-red.csv"))

#Looking for Null or Missing values in data
which(is.na(wine))

wine['logqua'] = ifelse(wine$quality >= 6.5 ,1,0)

#Frequentist Analysis
mmm <- glm(wine$logqua ~ wine$fixed.acidity + wine$volatile.acidity +  wine$citric.acid + wine$residual.sugar + wine$chlorides + wine$free.sulfur.dioxide + wine$total.sulfur.dioxide + wine$density + wine$pH + wine$sulphates + wine$alcohol, data = wine, family = binomial(link = 'logit') )
mmm$coefficients
tsulr <- seq(from = min(wine$total.sulfur.dioxide), to = max(wine$total.sulfur.dioxide), by = 0.1 )
newmod <- mmm$coefficients[1] + mmm$coefficients[2]*mean(wine$fixed.acidity) + mmm$coefficients[3]*mean(wine$volatile.acidity) + mmm$coefficients[4]*mean(wine$citric.acid) + mmm$coefficients[5]*mean(wine$residual.sugar) + mmm$coefficients[6]*mean(wine$chlorides) + mmm$coefficients[7]*mean(wine$free.sulfur.dioxide) + mmm$coefficients[8]*tsulr + mmm$coefficients[9]*mean(wine$density) + mmm$coefficients[10]*mean(wine$pH) + mmm$coefficients[11]*mean(wine$sulphates)+ mmm$coefficients[12]*mean(wine$alcohol)
newmodp <- exp(newmod)/(1 + exp(newmod))
plot(tsulr,newmodp/max(newmodp), ylim = c(0,1), xlab = "Amount of Total Sulphur Dioxide", ylab = "Probability of a Good Wine")
abline(h = 0.5, lty = 10)

  
#Bayesian Analysis
# Log-posterior distribution
logpost <- function(beta,x,y)
{
  asd <- as.numeric(x %*% beta)
  loga <- asd - log(1+exp(asd))
  logb <- log(1-exp(loga))
  logl <- sum(loga[y==1]) + sum(logb[y==0])
  lprior <- sum(dnorm(beta,0,10,log=T))
  return(logl + lprior)
}

S <- 10^4
X1=cbind(rep(1,nrow(wine)),wine$fixed.acidity,wine$volatile.acidity,wine$citric.acid,wine$residual.sugar,wine$chlorides,wine$free.sulfur.dioxide,wine$total.sulfur.dioxide,wine$density,wine$sulphates,wine$pH,wine$alcohol)
X = X1[1,]
X <- rbind(X,X1)
y =wine$logqua[1]
y <- c(y,wine$logqua)

#First initialisation *0.2
d1 <- mmm$coefficients*0.2
asd_mat1 <- matrix(NA,nrow=S,ncol=ncol(X))
asd_mat1[1,] <- as.numeric(d1)
sigma <- solve(t(X) %*% X)
k <- ncol(asd_mat1)
acc <- 1
for(iter in 2:S)
{
  bt1 <- rmvnorm(1,asd_mat1[iter-1,],sigma)
  newpost=logpost(t(bt1),X,y)
  oldpost=logpost(matrix(asd_mat1[iter-1,],ncol=1),X,y)
  
  if(runif(1,0,1)>exp(newpost-oldpost)){
    asd_mat1[iter,]=asd_mat1[iter-1,]
  } else{
    asd_mat1[iter,]=bt1
    acc=acc+1
  }
  
  if(iter%%1000==0){print(c(iter,acc/iter))}}

#Second initialisation *0.7
d2 <- (mmm$coefficients*0.7)
asd_mat2 <- matrix(NA,nrow=S,ncol=ncol(X))
asd_mat2[1,] <- as.numeric(d2)
sigma <- solve(t(X) %*% X)
k <- ncol(asd_mat2)
acc <- 1
for(iter in 2:S)
{
  bt2 <- rmvnorm(1,asd_mat2[iter-1,],sigma)
  newpost=logpost(t(bt2),X,y)
  oldpost=logpost(matrix(asd_mat2[iter-1,],ncol=1),X,y)
  
  if(runif(1,0,1)>exp(newpost-oldpost)){
    asd_mat2[iter,]=asd_mat2[iter-1,]
  } else{
    asd_mat2[iter,]=bt2
    acc=acc+1
  }
  
  if(iter%%1000==0){print(c(iter,acc/iter))}}

  #Third initialisation *1.2
  d3 <- (mmm$coefficients*1.2)
  asd_mat3 <- matrix(NA,nrow=S,ncol=ncol(X))
  asd_mat3[1,] <- as.numeric(d3)
  sigma <- solve(t(X) %*% X)
  k <- ncol(asd_mat3)
  acc <- 1
  for(iter in 2:S)
  {
    bt3 <- rmvnorm(1,asd_mat3[iter-1,],sigma)
    newpost=logpost(t(bt3),X,y)
    oldpost=logpost(matrix(asd_mat3[iter-1,],ncol=1),X,y)
    
    if(runif(1,0,1)>exp(newpost-oldpost)){
      asd_mat3[iter,]=asd_mat3[iter-1,]
    } else{
      asd_mat3[iter,]=bt3
      acc=acc+1
    }
    
    if(iter%%1000==0){print(c(iter,acc/iter))}}
 
  #Fourth initialisation *1.7
  d4 <- (mmm$coefficients*1.7)
  asd_mat4 <- matrix(NA,nrow=S,ncol=ncol(X))
  asd_mat4[1,] <- as.numeric(d4)
  sigma <- solve(t(X) %*% X)
  k <- ncol(asd_mat4)
  acc <- 1
  for(iter in 2:S)
  {
    bt4 <- rmvnorm(1,asd_mat4[iter-1,],sigma)
    newpost=logpost(t(bt4),X,y)
    oldpost=logpost(matrix(asd_mat4[iter-1,],ncol=1),X,y)
    
    if(runif(1,0,1)>exp(newpost-oldpost)){
      asd_mat4[iter,]=asd_mat4[iter-1,]
    } else{
      asd_mat4[iter,]=bt4
      acc=acc+1
    }

    if(iter%%1000==0){print(c(iter,acc/iter))}}
  
#Intercept Plots
  dev.new(width=5, height=4, unit="in")
  df1 <- data.frame(simulation = seq(1,10000, by =1), one = asd_mat1[,1],two = asd_mat2[,1],three = asd_mat3[,1], four = asd_mat4[,1],GLM = rep(mmm$coefficients[1],10000))
  ggplot(df1, aes(x=simulation))+ 
    geom_line(aes(y = one,color = "green")) +  
    geom_line(aes(y = two,color="blue")) +  
    geom_line(aes(y = three,color="red")) +    
    geom_line(aes(y = four,color="black"))+
    geom_line(aes(y = GLM,color="orange"))+
    scale_color_identity(name = "Initialisation", breaks = c("green", "blue", "red", "black", "orange"),
                         labels = c("1st", "2nd", "3rd", "4th", "GLM"),guide = "legend")+
    labs(y=expression(beta[0]),title = 'Intercept')
  
dev.new(width=5, height=4, unit="in")
par(mfrow=c(2,2))
plot(asd_mat1[,1],type="l", ylab=expression(beta[0]), main = 'Intercept1')
abline(h=d1[1],col="red",lty=2,lwd=5)
plot(asd_mat2[,1],type="l", ylab=expression(beta[0]), main = 'Intercept2')
abline(h=d2[1],col="red",lty=2,lwd=5)
plot(asd_mat3[,1],type="l", ylab=expression(beta[0]), main = 'Intercept3')
abline(h=d3[1],col="red",lty=2,lwd=5)
plot(asd_mat4[,1],type="l", ylab=expression(beta[0]), main = 'Intercept4')
abline(h=d4[1],col="red",lty=2,lwd=5)

#Fixed Acidity Plots
dev.new(width=5, height=4, unit="in")
df2 <- data.frame(simulation = seq(1,10000, by =1), one = asd_mat1[,2],two = asd_mat2[,2],three = asd_mat3[,2], four = asd_mat4[,2], GLM = rep(mmm$coefficients[2],10000))
ggplot(df2, aes(x=simulation))+ 
  geom_line(aes(y = one,color = "green")) +  
  geom_line(aes(y = two,color="blue")) +  
  geom_line(aes(y = three,color="red")) +    
  geom_line(aes(y = four,color="black"))+
  geom_line(aes(y = GLM,color="orange"))+
  scale_color_identity(name = "Initialisation", breaks = c("green", "blue", "red", "black", "orange"),
                       labels = c("1st", "2nd", "3rd", "4th", "GLM"),guide = "legend")+
  labs(y=expression(beta[1]),title = 'Fixed Acidity')

dev.new(width=5, height=4, unit="in")
par(mfrow=c(2,2))
plot(asd_mat1[,2],type="l", ylab=expression(beta[1]), main = 'Fixed Acidity 1')
abline(h=d1[2],col="red",lty=2,lwd=5)
plot(asd_mat2[,2],type="l", ylab=expression(beta[1]), main = 'Fixed Acidity 2')
abline(h=d2[2],col="red",lty=2,lwd=5)
plot(asd_mat3[,2],type="l", ylab=expression(beta[1]), main = 'Fixed Acidity 3')
abline(h=d3[2],col="red",lty=2,lwd=5)
plot(asd_mat4[,2],type="l", ylab=expression(beta[1]), main = 'Fixed Acidity 4')
abline(h=d4[2],col="red",lty=2,lwd=5)

#Volatile Acidity Plots
dev.new(width=5, height=4, unit="in")
df3 <- data.frame(simulation = seq(1,10000, by =1), one = asd_mat1[,3],two = asd_mat2[,3],three = asd_mat3[,3], four = asd_mat4[,3], GLM = rep(mmm$coefficients[3],10000))
ggplot(df3, aes(x=simulation))+ 
  geom_line(aes(y = one,color = "green")) +  
  geom_line(aes(y = two,color="blue")) +  
  geom_line(aes(y = three,color="red")) +    
  geom_line(aes(y = four,color="black"))+
  geom_line(aes(y = GLM,color="orange"))+
  scale_color_identity(name = "Initialisation", breaks = c("green", "blue", "red", "black", "orange"),
                       labels = c("1st", "2nd", "3rd", "4th", "GLM"),guide = "legend")+
  labs(y=expression(beta[2]),title = 'Volatile Acidity')

dev.new(width=5, height=4, unit="in")
par(mfrow=c(2,2))
plot(asd_mat1[,3],type="l", ylab=expression(beta[2]), main = 'Volatile Acidity 1')
abline(h=d1[3],col="red",lty=2,lwd=5)
plot(asd_mat2[,3],type="l", ylab=expression(beta[2]), main = 'Volatile Acidity 2')
abline(h=d2[3],col="red",lty=2,lwd=5)
plot(asd_mat3[,3],type="l", ylab=expression(beta[2]), main = 'Volatile Acidity 3')
abline(h=d3[3],col="red",lty=2,lwd=5)
plot(asd_mat4[,3],type="l", ylab=expression(beta[2]), main = 'Volatile Acidity 4')
abline(h=d4[3],col="red",lty=2,lwd=5)

#Citric Acid Plots
dev.new(width=5, height=4, unit="in")
df4 <- data.frame(simulation = seq(1,10000, by =1), one = asd_mat1[,4],two = asd_mat2[,4],three = asd_mat3[,4], four = asd_mat4[,4], GLM = rep(mmm$coefficients[4],10000))
ggplot(df4, aes(x=simulation))+ 
  geom_line(aes(y = one,color = "green")) +  
  geom_line(aes(y = two,color="blue")) +  
  geom_line(aes(y = three,color="red")) +    
  geom_line(aes(y = four,color="black"))+
  geom_line(aes(y = GLM,color="orange"))+
  scale_color_identity(name = "Initialisation", breaks = c("green", "blue", "red", "black", "orange"),
                       labels = c("1st", "2nd", "3rd", "4th", "GLM"),guide = "legend")+
  labs(y=expression(beta[3]),title = 'Citric Acid')

dev.new(width=5, height=4, unit="in")
par(mfrow=c(2,2))
plot(asd_mat1[,4],type="l", ylab=expression(beta[3]), main = 'Citric Acid 1')
abline(h=d1[4],col="red",lty=2,lwd=5)
plot(asd_mat2[,4],type="l", ylab=expression(beta[3]), main = 'Citric Acid 2')
abline(h=d2[4],col="red",lty=2,lwd=5)
plot(asd_mat3[,4],type="l", ylab=expression(beta[3]), main = 'Citric Acid 3')
abline(h=d3[4],col="red",lty=2,lwd=5)
plot(asd_mat4[,4],type="l", ylab=expression(beta[3]), main = 'Citric Acid 4')
abline(h=d4[4],col="red",lty=2,lwd=5)

#Residual Sugar
dev.new(width=5, height=4, unit="in")
df5 <- data.frame(simulation = seq(1,10000, by =1), one = asd_mat1[,5],two = asd_mat2[,5],three = asd_mat3[,5], four = asd_mat4[,5], GLM = rep(mmm$coefficients[5],10000))
ggplot(df5, aes(x=simulation))+ 
  geom_line(aes(y = one,color = "green")) +  
  geom_line(aes(y = two,color="blue")) +  
  geom_line(aes(y = three,color="red")) +    
  geom_line(aes(y = four,color="black"))+
  geom_line(aes(y = GLM,color="orange"))+
  scale_color_identity(name = "Initialisation", breaks = c("green", "blue", "red", "black", "orange"),
                       labels = c("1st", "2nd", "3rd", "4th", "GLM"),guide = "legend")+
  labs(y=expression(beta[4]),title = 'Residual Sugar')

dev.new(width=5, height=4, unit="in")
par(mfrow=c(2,2))
plot(asd_mat1[,5],type="l", ylab=expression(beta[4]), main = 'Residual Sugar 1')
abline(h=d1[5],col="red",lty=2,lwd=5)
plot(asd_mat2[,5],type="l", ylab=expression(beta[4]), main = 'Residual Sugar 2')
abline(h=d2[5],col="red",lty=2,lwd=5)
plot(asd_mat3[,5],type="l", ylab=expression(beta[4]), main = 'Residual Sugar 3')
abline(h=d3[5],col="red",lty=2,lwd=5)
plot(asd_mat4[,5],type="l", ylab=expression(beta[4]), main = 'Residual Sugar 4')
abline(h=d4[5],col="red",lty=2,lwd=5)

#Chlorides
dev.new(width=5, height=4, unit="in")
df6 <- data.frame(simulation = seq(1,10000, by =1), one = asd_mat1[,6],two = asd_mat2[,6],three = asd_mat3[,6], four = asd_mat4[,6], GLM = rep(mmm$coefficients[6],10000))
ggplot(df6, aes(x=simulation))+ 
  geom_line(aes(y = one,color = "green")) +  
  geom_line(aes(y = two,color="blue")) +  
  geom_line(aes(y = three,color="red")) +    
  geom_line(aes(y = four,color="black"))+
  geom_line(aes(y = GLM,color="orange"))+
  scale_color_identity(name = "Initialisation", breaks = c("green", "blue", "red", "black", "orange"),
                       labels = c("1st", "2nd", "3rd", "4th", "GLM"),guide = "legend")+
  labs(y=expression(beta[5]),title = 'Chlorides')

dev.new(width=5, height=4, unit="in")
par(mfrow=c(2,2))
plot(asd_mat1[,6],type="l", ylab=expression(beta[5]), main = 'Chlorides 1')
abline(h=d1[6],col="red",lty=2,lwd=5)
plot(asd_mat2[,6],type="l", ylab=expression(beta[5]), main = 'Chlorides 2')
abline(h=d2[6],col="red",lty=2,lwd=5)
plot(asd_mat3[,6],type="l", ylab=expression(beta[5]), main = 'Chlorides 3')
abline(h=d3[6],col="red",lty=2,lwd=5)
plot(asd_mat4[,6],type="l", ylab=expression(beta[5]), main = 'Chlorides 4')
abline(h=d4[6],col="red",lty=2,lwd=5)

#Free Sulfur Dioxide
dev.new(width=5, height=4, unit="in")
df7 <- data.frame(simulation = seq(1,10000, by =1), one = asd_mat1[,7],two = asd_mat2[,7],three = asd_mat3[,7], four = asd_mat4[,7], GLM = rep(mmm$coefficients[7],10000))
ggplot(df7, aes(x=simulation))+ 
  geom_line(aes(y = one,color = "green")) +  
  geom_line(aes(y = two,color="blue")) +  
  geom_line(aes(y = three,color="red")) +    
  geom_line(aes(y = four,color="black"))+
  geom_line(aes(y = GLM,color="orange"))+
  scale_color_identity(name = "Initialisation", breaks = c("green", "blue", "red", "black", "orange"),
                       labels = c("1st", "2nd", "3rd", "4th", "GLM"),guide = "legend")+
  labs(y=expression(beta[6]),title = 'Free Sulfur Dioxide')

dev.new(width=5, height=4, unit="in")
par(mfrow=c(2,2))
plot(asd_mat1[,7],type="l", ylab=expression(beta[6]), main = 'Free Sulfur Dioxide 1')
abline(h=d1[7],col="red",lty=2,lwd=5)
plot(asd_mat2[,7],type="l", ylab=expression(beta[6]), main = 'Free Sulfur Dioxide 2')
abline(h=d2[7],col="red",lty=2,lwd=5)
plot(asd_mat3[,7],type="l", ylab=expression(beta[6]), main = 'Free Sulfur Dioxide 3')
abline(h=d3[7],col="red",lty=2,lwd=5)
plot(asd_mat4[,7],type="l", ylab=expression(beta[0]), main = 'Free Sulfur Dioxide 4')
abline(h=d4[7],col="red",lty=2,lwd=5)

#Total Sulfur Dioxide
dev.new(width=5, height=4, unit="in")
df8 <- data.frame(simulation = seq(1,10000, by =1), one = asd_mat1[,8],two = asd_mat2[,8],three = asd_mat3[,8], four = asd_mat4[,8], GLM = rep(mmm$coefficients[8],10000))
ggplot(df8, aes(x=simulation))+ 
  geom_line(aes(y = one,color = "green")) +  
  geom_line(aes(y = two,color="blue")) +  
  geom_line(aes(y = three,color="red")) +    
  geom_line(aes(y = four,color="black"))+
  geom_line(aes(y = GLM,color="orange"))+
  scale_color_identity(name = "Initialisation", breaks = c("green", "blue", "red", "black", "orange"),
                       labels = c("1st", "2nd", "3rd", "4th", "GLM"),guide = "legend")+
  labs(y=expression(beta[7]),title = 'Total Sulfur Dioxide')

dev.new(width=5, height=4, unit="in")
par(mfrow=c(2,2))
plot(asd_mat1[,8],type="l", ylab=expression(beta[7]), main = 'Total Sulfur Dioxide 1')
abline(h=d1[8],col="red",lty=2,lwd=5)
plot(asd_mat2[,8],type="l", ylab=expression(beta[7]), main = 'Total Sulfur Dioxide 2')
abline(h=d2[8],col="red",lty=2,lwd=5)
plot(asd_mat3[,8],type="l", ylab=expression(beta[7]), main = 'Total Sulfur Dioxide 3')
abline(h=d3[8],col="red",lty=2,lwd=5)
plot(asd_mat4[,8],type="l", ylab=expression(beta[7]), main = 'Total Sulfur Dioxide 4')
abline(h=d4[8],col="red",lty=2,lwd=5)

#Density
dev.new(width=5, height=4, unit="in")
df9 <- data.frame(simulation = seq(1,10000, by =1), one = asd_mat1[,9],two = asd_mat2[,9],three = asd_mat3[,9], four = asd_mat4[,9], GLM = rep(mmm$coefficients[9],10000))
ggplot(df9, aes(x=simulation))+ 
  geom_line(aes(y = one,color = "green")) +  
  geom_line(aes(y = two,color="blue")) +  
  geom_line(aes(y = three,color="red")) +    
  geom_line(aes(y = four,color="black"))+
  geom_line(aes(y = GLM,color="orange"))+
  scale_color_identity(name = "Initialisation", breaks = c("green", "blue", "red", "black", "orange"),
                       labels = c("1st", "2nd", "3rd", "4th", "GLM"),guide = "legend")+
  labs(y=expression(beta[8]),title = 'Density')

dev.new(width=5, height=4, unit="in")
par(mfrow=c(2,2))
plot(asd_mat1[,9],type="l", ylab=expression(beta[8]), main = 'Density 1')
abline(h=d1[9],col="red",lty=2,lwd=5)
plot(asd_mat2[,9],type="l", ylab=expression(beta[8]), main = 'Density 2')
abline(h=d2[9],col="red",lty=2,lwd=5)
plot(asd_mat3[,9],type="l", ylab=expression(beta[8]), main = 'Density 3')
abline(h=d3[9],col="red",lty=2,lwd=5)
plot(asd_mat4[,9],type="l", ylab=expression(beta[8]), main = 'Density 4')
abline(h=d4[9],col="red",lty=2,lwd=5)

#pH 
dev.new(width=5, height=4, unit="in")
df10 <- data.frame(simulation = seq(1,10000, by =1), one = asd_mat1[,10],two = asd_mat2[,10],three = asd_mat3[,10], four = asd_mat4[,10], GLM = rep(mmm$coefficients[10],10000))
ggplot(df10, aes(x=simulation))+ 
  geom_line(aes(y = one,color = "green")) +  
  geom_line(aes(y = two,color="blue")) +  
  geom_line(aes(y = three,color="red")) +    
  geom_line(aes(y = four,color="black"))+
  geom_line(aes(y = GLM,color="orange"))+
  scale_color_identity(name = "Initialisation", breaks = c("green", "blue", "red", "black", "orange"),
                       labels = c("1st", "2nd", "3rd", "4th", "GLM"),guide = "legend")+
  labs(y=expression(beta[9]),title = 'pH')

dev.new(width=5, height=4, unit="in")
par(mfrow=c(2,2))
plot(asd_mat1[,10],type="l", ylab=expression(beta[9]), main = 'pH 1')
abline(h=d1[10],col="red",lty=2,lwd=5)
plot(asd_mat2[,10],type="l", ylab=expression(beta[9]), main = 'pH 2')
abline(h=d2[10],col="red",lty=2,lwd=5)
plot(asd_mat3[,10],type="l", ylab=expression(beta[9]), main = 'pH 3')
abline(h=d3[10],col="red",lty=2,lwd=5)
plot(asd_mat4[,10],type="l", ylab=expression(beta[9]), main = 'pH 4')
abline(h=d4[10],col="red",lty=2,lwd=5)

#Sulphates 
dev.new(width=5, height=4, unit="in")
df11 <- data.frame(simulation = seq(1,10000, by =1), one = asd_mat1[,11],two = asd_mat2[,11],three = asd_mat3[,11], four = asd_mat4[,11], GLM = rep(mmm$coefficients[11],10000))
ggplot(df11, aes(x=simulation))+ 
  geom_line(aes(y = one,color = "green")) +  
  geom_line(aes(y = two,color="blue")) +  
  geom_line(aes(y = three,color="red")) +    
  geom_line(aes(y = four,color="black"))+
  geom_line(aes(y = GLM,color="orange"))+
  scale_color_identity(name = "Initialisation", breaks = c("green", "blue", "red", "black", "orange"),
                       labels = c("1st", "2nd", "3rd", "4th", "GLM"),guide = "legend")+
  labs(y=expression(beta[10]),title = 'Sulphates')

dev.new(width=5, height=4, unit="in")
par(mfrow=c(2,2))
plot(asd_mat1[,11],type="l", ylab=expression(beta[10]), main = 'Sulphates 1')
abline(h=d1[11],col="red",lty=2,lwd=5)
plot(asd_mat2[,11],type="l", ylab=expression(beta[10]), main = 'Sulphates 2')
abline(h=d2[11],col="red",lty=2,lwd=5)
plot(asd_mat3[,11],type="l", ylab=expression(beta[10]), main = 'Sulphates 3')
abline(h=d3[11],col="red",lty=2,lwd=5)
plot(asd_mat4[,11],type="l", ylab=expression(beta[10]), main = 'Sulphates 4')
abline(h=d4[11],col="red",lty=2,lwd=5)

#Alcohol 
dev.new(width=5, height=4, unit="in")
df12 <- data.frame(simulation = seq(1,10000, by =1), one = asd_mat1[,12],two = asd_mat2[,12],three = asd_mat3[,12], four = asd_mat4[,12], GLM = rep(mmm$coefficients[12],10000))
ggplot(df12, aes(x=simulation))+ 
  geom_line(aes(y = one,color = "green")) +  
  geom_line(aes(y = two,color="blue")) +  
  geom_line(aes(y = three,color="red")) +    
  geom_line(aes(y = four,color="black"))+
  geom_line(aes(y = GLM,color="orange"))+
  scale_color_identity(name = "Initialisation", breaks = c("green", "blue", "red", "black", "orange"),
                       labels = c("1st", "2nd", "3rd", "4th", "GLM"),guide = "legend")+
  labs(y=expression(beta[11]),title = 'Alcohol')

dev.new(width=5, height=4, unit="in")
par(mfrow=c(2,2))
plot(asd_mat1[,12],type="l", ylab=expression(beta[11]), main = 'Alcohol 1')
abline(h=d1[12],col="red",lty=2,lwd=5)
plot(asd_mat2[,12],type="l", ylab=expression(beta[11]), main = 'Alcohol 2')
abline(h=d2[12],col="red",lty=2,lwd=5)
plot(asd_mat3[,12],type="l", ylab=expression(beta[11]), main = 'Alcohol 3')
abline(h=d3[12],col="red",lty=2,lwd=5)
plot(asd_mat4[,12],type="l", ylab=expression(beta[11]), main = 'Alcohol 4')
abline(h=d4[12],col="red",lty=2,lwd=5)

pre<- c(asd_mat1[iter,], asd_mat2[iter,], asd_mat3[iter,], asd_mat4[iter,])

#Prediction 
x_new <-c(7.5,0.6,0,1.7,0.085,5,45,0.9965,3.40,0.63,12)
p_new1 <- exp(mmm$coefficients[2:12] * x_new) / (1 + exp(mmm$coefficients[2:12] * x_new))
p_new2 <- exp(pre[2:12] * x_new) / (1 + exp(pre[2:12] * x_new))

dfn <-data.frame(Covariate =rep(c('Fixed Acidity ','Volatile Acidity','Citric Acid','Residual Sugar','Chlorides','Free SO2','Total SO2','Density','pH','Sulphates','Alcohol'),2), model = rep(c('Frequentist', 'Bayesian'),each =11),'Probability of a Good Wine' = c(p_new1, p_new2))
dev.new(width=5, height=4, unit="in")
ggplot(dfn, aes(x=Covariate, y=Probability.of.a.Good.Wine, fill=model)) + 
geom_bar(stat="identity", position=position_dodge())

#Metrop Analysis
lat <- function(x, y) function(beta)
{
  asd <- as.numeric(x %*% beta)
  loga <- asd - log(1+exp(asd))
  logb <- log(1-exp(loga))
  logl <- sum(loga[y==1]) + sum(logb[y==0])
  lprior <- sum(dnorm(beta,0,10,log=T))
  return(logl + lprior)
}
lut <- lat(X1, mmm$y)
coff <- as.numeric(coefficients(mmm))
zxc <- metrop(lut, coff, S,blen = 1,nspac = 10,debug = TRUE)

#Comparison Plots
#Intercept
dev.new(width=5, height=4, unit="in")
dfm1 <- data.frame(simulation = seq(1,10000, by =1), mha = asd_mat1[,1],met = zxc$batch[,1],GLM = rep(mmm$coefficients[1],10000))
ggplot(dfm1, aes(x=simulation)) + 
  geom_line(aes(y = mha,color = 'black')) + 
  geom_line(aes(y = met,color="navyblue"), lwd =1.5)+
  geom_line(aes(y = GLM,color="red"),lwd =1.5)+ 
  scale_color_identity(name = "Model", breaks = c("black", "navyblue", "red"),
  labels = c("Bayesian", "Metrop", "Frequentist"),guide = "legend")+
  labs(y=expression(beta[0]),title = 'Intercept')

#Fixed Acidity
dev.new(width=5, height=4, unit="in")
dfm2 <- data.frame(simulation = seq(1,10000, by =1), mha = asd_mat1[,2],met = zxc$batch[,2],GLM = rep(mmm$coefficients[2],10000))
ggplot(dfm2, aes(x=simulation)) + 
  geom_line(aes(y = mha,color = 'black')) + 
  geom_line(aes(y = met,color="navyblue"), lwd =1.5)+
  geom_line(aes(y = GLM,color="red"),lwd =1.5)+ 
  scale_color_identity(name = "Model", breaks = c("black", "navyblue", "red"),
  labels = c("Bayesian", "Metrop", "Frequentist"),guide = "legend")+
  labs(y=expression(beta[0]),title = 'Fixed Acidity')

#Volatile Acidity
dev.new(width=5, height=4, unit="in")
dfm3 <- data.frame(simulation = seq(1,10000, by =1), mha = asd_mat1[,3],met = zxc$batch[,3],GLM = rep(mmm$coefficients[3],10000))
ggplot(dfm3, aes(x=simulation)) + 
  geom_line(aes(y = mha,color = 'black')) + 
  geom_line(aes(y = met,color="navyblue"), lwd =1.5)+
  geom_line(aes(y = GLM,color="red"),lwd =1.5)+ 
  scale_color_identity(name = "Model", breaks = c("black", "navyblue", "red"),
  labels = c("Bayesian", "Metrop", "Frequentist"),guide = "legend")+
  labs(y=expression(beta[0]),title = 'Volatile Acidity')

#Citric Acid
dev.new(width=5, height=4, unit="in")
dfm4 <- data.frame(simulation = seq(1,10000, by =1), mha = asd_mat1[,4],met = zxc$batch[,4],GLM = rep(mmm$coefficients[4],10000))
ggplot(dfm4, aes(x=simulation)) + 
  geom_line(aes(y = mha,color = 'black')) + 
  geom_line(aes(y = met,color="navyblue"), lwd =1.5)+
  geom_line(aes(y = GLM,color="red"),lwd =1.5)+ 
  scale_color_identity(name = "Model", breaks = c("black", "navyblue", "red"),
  labels = c("Bayesian", "Metrop", "Frequentist"),guide = "legend")+
  labs(y=expression(beta[0]),title = 'Citric Acid')

#Residual Sugar
dev.new(width=5, height=4, unit="in")
dfm5 <- data.frame(simulation = seq(1,10000, by =1), mha = asd_mat1[,5],met = zxc$batch[,5],GLM = rep(mmm$coefficients[5],10000))
ggplot(dfm5, aes(x=simulation)) + 
  geom_line(aes(y = mha,color = 'black')) + 
  geom_line(aes(y = met,color="navyblue"), lwd =1.5)+
  geom_line(aes(y = GLM,color="red"),lwd =1.5)+ 
  scale_color_identity(name = "Model", breaks = c("black", "navyblue", "red"),
  labels = c("Bayesian", "Metrop", "Frequentist"),guide = "legend")+
  labs(y=expression(beta[0]),title = 'Residual Sugar')

#Chlorides
dev.new(width=5, height=4, unit="in")
dfm6 <- data.frame(simulation = seq(1,10000, by =1), mha = asd_mat1[,6],met = zxc$batch[,6],GLM = rep(mmm$coefficients[6],10000))
ggplot(dfm6, aes(x=simulation)) + 
  geom_line(aes(y = mha,color = 'black')) + 
  geom_line(aes(y = met,color="navyblue"), lwd =1.5)+
  geom_line(aes(y = GLM,color="red"),lwd =1.5)+ 
  scale_color_identity(name = "Model", breaks = c("black", "navyblue", "red"),
  labels = c("Bayesian", "Metrop", "Frequentist"),guide = "legend")+
  labs(y=expression(beta[0]),title = 'Chlorides')

#Free Sulfur Dioxide
dev.new(width=5, height=4, unit="in")
dfm7 <- data.frame(simulation = seq(1,10000, by =1), mha = asd_mat1[,7],met = zxc$batch[,7],GLM = rep(mmm$coefficients[7],10000))
ggplot(dfm7, aes(x=simulation)) + 
  geom_line(aes(y = mha,color = 'black')) + 
  geom_line(aes(y = met,color="navyblue"), lwd =1.5)+
  geom_line(aes(y = GLM,color="red"),lwd =1.5)+ 
  scale_color_identity(name = "Model", breaks = c("black", "navyblue", "red"),
  labels = c("Bayesian", "Metrop", "Frequentist"),guide = "legend")+
  labs(y=expression(beta[0]),title = 'Free Sulfur Dioxide')

#Total Sulfur Dioxide
dev.new(width=5, height=4, unit="in")
dfm8 <- data.frame(simulation = seq(1,10000, by =1), mha = asd_mat1[,8],met = zxc$batch[,8],GLM = rep(mmm$coefficients[8],10000))
ggplot(dfm8, aes(x=simulation)) + 
  geom_line(aes(y = mha,color = 'black')) + 
  geom_line(aes(y = met,color="navyblue"), lwd =1.5)+
  geom_line(aes(y = GLM,color="red"),lwd =1.5)+ 
  scale_color_identity(name = "Model", breaks = c("black", "navyblue", "red"),
  labels = c("Bayesian", "Metrop", "Frequentist"),guide = "legend")+
  labs(y=expression(beta[0]),title = 'Total Sulfur Dioxide')

#Density
dev.new(width=5, height=4, unit="in")
dfm9 <- data.frame(simulation = seq(1,10000, by =1), mha = asd_mat1[,9],met = zxc$batch[,9],GLM = rep(mmm$coefficients[9],10000))
ggplot(dfm9, aes(x=simulation)) + 
  geom_line(aes(y = mha,color = 'black')) + 
  geom_line(aes(y = met,color="navyblue"), lwd =1.5)+
  geom_line(aes(y = GLM,color="red"),lwd =1.5)+ 
  scale_color_identity(name = "Model", breaks = c("black", "navyblue", "red"),
  labels = c("Bayesian", "Metrop", "Frequentist"),guide = "legend")+
  labs(y=expression(beta[0]),title = 'Density')

#pH
dev.new(width=5, height=4, unit="in")
dfm10 <- data.frame(simulation = seq(1,10000, by =1), mha = asd_mat1[,10],met = zxc$batch[,10],GLM = rep(mmm$coefficients[10],10000))
ggplot(dfm10, aes(x=simulation)) + 
  geom_line(aes(y = mha,color = 'black')) + 
  geom_line(aes(y = met,color="navyblue"), lwd =1.5)+
  geom_line(aes(y = GLM,color="red"),lwd =1.5)+ 
  scale_color_identity(name = "Model", breaks = c("black", "navyblue", "red"),
  labels = c("Bayesian", "Metrop", "Frequentist"),guide = "legend")+
  labs(y=expression(beta[0]),title = 'pH')

#Sulphates
dev.new(width=5, height=4, unit="in")
dfm11 <- data.frame(simulation = seq(1,10000, by =1), mha = asd_mat1[,11],met = zxc$batch[,11],GLM = rep(mmm$coefficients[11],10000))
ggplot(dfm11, aes(x=simulation)) + 
  geom_line(aes(y = mha,color = 'black')) + 
  geom_line(aes(y = met,color="navyblue"), lwd =1.5)+
  geom_line(aes(y = GLM,color="red"),lwd =1.5)+ 
  scale_color_identity(name = "Model", breaks = c("black", "navyblue", "red"),
  labels = c("Bayesian", "Metrop", "Frequentist"),guide = "legend")+
  labs(y=expression(beta[0]),title = 'Sulphates')

#Alcohol
dev.new(width=5, height=4, unit="in")
dfm12 <- data.frame(simulation = seq(1,10000, by =1), mha = asd_mat1[,12],met = zxc$batch[,12],GLM = rep(mmm$coefficients[12],10000))
ggplot(dfm12, aes(x=simulation)) + 
  geom_line(aes(y = mha,color = 'black')) + 
  geom_line(aes(y = met,color="navyblue"), lwd =1.5)+
  geom_line(aes(y = GLM,color="red"),lwd =1.5)+ 
  scale_color_identity(name = "Model", breaks = c("black", "navyblue", "red"),
  labels = c("Bayesian", "Metrop", "Frequentist"),guide = "legend")+
  labs(y=expression(beta[0]),title = 'Alcohol')