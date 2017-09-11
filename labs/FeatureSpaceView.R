
nfeat = 1000
lambda = 0.1

phi <- function(x){
  cc <- seq(0,1,length=nfeat)
  exp(-(x-cc)^2/(2*lambda^2))
}



x <- seq(0,1,len=7)
y = x*sin(x*5)
plot(x,y)
pdnorm <-function(x, cc){
  0.3*exp(-(x-cc)^2/(2*lambda^2))-1
}

xp <- seq(0,1,0.001)
for(cc in seq(0,1,len=min(10, nfeat))){
  lines(xp, pdnorm(xp, cc), lwd=0.3)
}


X <-t(sapply(x, phi))
sigma2 <- 0.00001
betahat <- solve(t(X)%*%X+sigma2*diag(dim(X)[2]) )%*%t(X)%*%(matrix(y, nc=1))

yp = t(sapply(xp, phi))%*%betahat
lines(xp, yp)

