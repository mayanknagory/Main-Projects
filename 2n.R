library(dplyr)
#TARGET DISTRIBUTION

set.seed(1234)
c<-seq(0,1, length.out =10000)
rtarget <- (0.33*rbeta(c,0,5)+0.33*rbeta(c,3,5)+0.33*rbeta(c,11,5))
dtarget <- (0.33*dbeta(c,0,5)+0.33*dbeta(c,3,5)+0.33*dbeta(c,11,5))

#PLOTTING TARGET FUNCTION
dev.new(width=5, height=4, unit="in")
plot(c,dtarget,xlab ='x', ylab ='f(x)',main = "Target Function")

#MEAN AND MODE OF TARGET FUNCTION
df <- data.frame(x = c, f = 0.33*dbeta(c,0,5)+0.33*dbeta(c,3,5)+0.33*dbeta(c,11,5))
df<- df %>% arrange(desc(f))
df[1:3,]

#MODE OF TARGET FUNCTION 
K <- df[2,1]
K
# K = 0.6942694

#MEAN OF TARGET FUNCTION
tmax <- mean(rtarget)
tmax

#Theoretical Acceptance Rate
1/(0.33*dbeta(K,0,5)+0.33*dbeta(K,3,5)+0.33*dbeta(K,11,5))


#ACCEPTANCE AND REJECTION ALGORITHM
n = 10000
theta = runif(n)
asd=runif(n,0,0.33*dbeta(K,0,5)+0.33*dbeta(K,3,5)+0.33*dbeta(K,11,5))
qwe= asd<0.33*dbeta(theta,0,5)+0.33*dbeta(theta,3,5)+0.33*dbeta(theta,11,5)
sum(qwe ==T)/n
mean(theta[qwe])

#IMPORTANCE SAMPLING
w = 0.33*dbeta(K,0,5)+0.33*dbeta(K,3,5)+0.33*dbeta(K,11,5)/dunif(theta)
wmean = sum(w*theta)/sum(w)
zxc = w/sum(w)
wmean

#GRAPHS AND PLOTS
dev.new(width=5, height=4, unit="in")
par(mfrow=c(2,2))
  
plot(theta,asd,xlab='x',pch=20,cex=0.5,ylab="g(x)",ylim=c(0,1.4))
lines(c,dtarget,col="red",lwd=3)
lines(c(0,1),c(0.33*dbeta(K,0,5)+0.33*dbeta(K,3,5)+0.33*dbeta(K,11,5),0.33*dbeta(K,0,5)+0.33*dbeta(K,3,5)+0.33*dbeta(K,11,5)),lwd=3,col="green")

plot(theta[qwe],asd[qwe],xlab='x',pch=20,
     cex=0.5,ylab="g(x)",xlim=c(0,1),ylim=c(0,1.4))
lines(c,dtarget,col="red",lwd=3)
lines(c(0,1),c(0.33*dbeta(K,0,5)+0.33*dbeta(K,3,5)+0.33*dbeta(K,11,5),0.33*dbeta(K,0,5)+0.33*dbeta(K,3,5)+0.33*dbeta(K,11,5)),lwd=3,col="green")

hist(theta[qwe],probability=T,xlab='x',
     ylab="Density",main="",xlim=c(0,1),ylim=c(0,2))
lines(c,dtarget,col="red",lwd=3)

hist(theta[qwe],probability=T,xlab='x',ylim=c(0,2),ylab="Density",main="")
lines(c,dtarget,lwd=5,col="blue")
d=density(theta,weights=zxc,from=0,to=1)
lines(d,col="cyan",lwd=5,lty=2)  
legend("topleft", legend=c("Accept/Rejection", "Importance Sampling"),col=c("blue","cyan"),lty = 1:1)