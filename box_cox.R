##box cox trans.
##2017.12.2
## question unsolved :
##    1. uniroot interval: is (-3,3) enough? (maybe use log.mle and optim)
##    2. de.log.mle may have same sign in -3,3 ,causing uniroot failed (maybe use log.mle and optim)
##    3. not a efficient way to find hat matrix and mse of lm model
library(magrittr)


 

# argument: lm model
# return : lambda, transformed y, new lm model
bc.transform = function(lm){
  if(class(lm)!="lm") return("put in fucking lm model only")
  # get y and x matrix
  y = lm$model[,1]
  x = lm$model
  x[,1] =1
  x = as.matrix(x)
  #save y names
  y.names = names(lm$model)[1]
  
  
  #hat matrix
  hat.x = x %*% solve(t(x) %*% x) %*% t(x)
  i.minus.hat = diag(dim(hat.x)[1]) - hat.x
  
  #box cox transform: find y and de.y
  bc.trans = function(y,lambda){
    if(lambda!=0){
      bc = (y^lambda-1)/lambda
      de.bc = (lambda*(y^lambda)*log(y)-y^lambda+1)/lambda^2
      return(list(y=bc,de.y=de.bc))
    } 
    if(lambda==0){
      bc = log(y)
      de.bc = log(y)^2/2
      return(list(y=bc,de.y=de.bc))    
    }
  }
  #de.lambda.mle
  de.log.mle = function(y,lambda){
    bc=bc.trans(y,lambda)
    lm.temp = lm$model
    lm.temp[,1] = bc$y
    names(lm.temp)[1] ="y"
    mean.sq = lm(y~.,data = lm.temp)$residuals^2 %>% mean()
    #mean.sq = (bc$y %*%i.minus.hat %*%bc$y)/ (length(bc$y)- dim(lm.temp)[2]) %>% as.vector()
    #print(mean.sq)
    return(sum(log(y))-(bc$de.y %*% i.minus.hat %*% bc$y)/mean.sq)
  }
  #find uniroot
  lambda = uniroot(de.log.mle,y=y,c(-3,3),tol=0.0001)$root
  
  #fit new model
  trans.data =lm$model
  new.y =  bc.trans(y,lambda)$y
  trans.data[,1] = new.y
  names(trans.data)[1]="y"
  lm.new = lm(y~.,data = trans.data)
  names(lm.new$model)[1]= y.names
  #res
  res = list(lambda=lambda,y.transform=new.y,lm=lm.new)
  class(res) ="bc"
  return(res)
  
}

### s3 method
print.bc = function(bc){
  print(bc[1]);print(bc[3])
}

qqplot = function(bc){
  UseMethod("qqplot")
}

qqplot.bc = function(bc){
  rstd.res = rstandard(bc$lm)
  qqnorm(rstd.res)
  qqline(rstd.res,col="red")
}

plot.bc = function(bc){
  plot(bc$lm)
}

### test data
###

set.seed(689)
x1 = rnorm(200,50,13)
x2 = rnorm(200,15,3)
x3 = rnorm(200,60,7)
y = x1^2-x2-x3 - rpois(200,l = 30)
lm2 = lm(y~x1+x2-x3)
temp = MASS::boxcox(lm2,l=seq(0,1,0.01))
temp$x[which.max(temp$y)] 

bc.transform(lm2)$lambda



### for HW
df = read.csv("jj.csv")
lm1 = lm(BIO~SAL+pH+K+Na+Zn,data=df)
rstandard(lm1) %T>% qqnorm() %>%  qqline(.,col="red")

temp = MASS::boxcox(BIO~SAL+pH+K+Na+Zn,data=df,l=seq(-2,2,0.001)) 
temp$x[which.max(temp$y)]

#
res = bc.transform(lm1)









