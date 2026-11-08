### Fermi-Dirac Probability Distribution

# Mehdi Shams
# mehdishams@kashanu.ac.ir

# Mohammad Ali Mirzaie
# m.a.mirzaie@ut.ac.ir

# In this session, we write some functions for calculating
# the measurements based on the theory of the papers


##### densities #####

### density of FD type I
density.FD1=function(x,beta,a=0){
  (exp(-beta*x))/
    (((exp(-beta*a))/(exp(beta)-1))+exp(-beta*x))
}
### density of FD type II
density.FD2=function(x,beta,a=0){
  (beta*exp(-beta*x))/
    (log(1+exp(-beta*a),base=exp(1))*(1+exp(-beta*x)))
}

# In both functions, x is the input of the density function,
# beta is the parameter, and a is the location parameter

# example
density.FD1(x = 1, beta = 1, a = 0)
density.FD2(x = 1, beta = 1, a = 0)

##### moments type R integral #####
### FD type one
moments.FD1_1=function(beta,a,k){
  f=function(x,beta,a,k){
  (x^k)*(exp(-beta*x))/
      (((exp(-beta*a))/(exp(beta)-1))+exp(-beta*x))
  }
  integrate(f,beta,a,k,lower=a,upper=Inf)
}
### FD type two
moments.FD2_1=function(beta,a,k){
  f=function(x,beta,a,k){
  (x^k)*(beta*exp(-beta*x))/
      (log(1+exp(-beta*a),base=exp(1))*(1+exp(-beta*x)))
  }
  integrate(f,beta,a,k,lower=a,upper=Inf)
}

# beta is the parameter, a is the location parameter,
# and k is the order of the moments.

# example
moments.FD1_1(beta = 0.25, a = 0, k = 1)
moments.FD1_1(beta = 0.25, a = 0, k = 2)

moments.FD2_1(beta = 0.25, a = 0, k = 1)
moments.FD2_1(beta = 0.25, a = 0, k = 2)


##### moments type MCMC #####
### FD type one
moments.FD1_2=function(beta,a,k,n){
  x=rexp(n,rate=beta)
  x=x+a
  y=(x^k)*(1/beta)*exp(-beta*a)/
    (((exp(-beta*a))/(exp(beta)-1))+exp(-beta*x))
  mean(y)
}
### FD type two
moments.FD2_2=function(beta,a,k,n){
  x=rexp(n,rate=beta)
  x=x+a
  y=(x^k)*exp(-beta*a)/
    (log(1+exp(-beta*a),base=exp(1))*(1+exp(-beta*x)))
  mean(y)
}

# beta is the parameter, a is the location parameter, k is
# the order of the moments, and n is the size of sample for MCMC

# example
moments.FD1_2(beta = 0.25, a = 0, k = 1, n = 1e6)
moments.FD1_2(beta = 0.25, a = 0, k = 2, n = 1e6)

moments.FD2_2(beta = 0.25, a = 0, k = 1, n = 1e6)
moments.FD2_2(beta = 0.25, a = 0, k = 2, n = 1e6)

##### fisher information R integral #####
### FD type one
fishinf.FD1_1=function(beta){
  f=function(x,beta){
  ((-x-(((-exp(beta)/(exp(beta)-1)^2 )-x*exp(-beta*x))/
          ((1/(exp(beta)-1))+exp(-beta*x)) ))^2)*
      (exp(-beta*x)/((1/(exp(beta)-1))+exp(-beta*x)))
  }
  integrate(f,beta,lower=0,upper=Inf)
}
### FD type two
fishinf.FD2_1=function(beta){
  f=function(x,beta){
  (((1/beta)-x+(x*exp(-beta*x))/(exp(-beta*x)+1))^2)*
      (beta*exp(-beta*x))/((exp(-beta*x)+1)*log(2,base=exp(1)))
  }
  integrate(f,beta,lower=0,upper=Inf)
}

#example
fishinf.FD1_1(beta = 0.25)
fishinf.FD2_1(beta = 0.25)

##### fisher information MCMC #####
### FD type one
fishinf.FD1_2=function(beta,n){
  x=rexp(n,rate=beta)
  y=((-x-(((-exp(beta)/(exp(beta)-1)^2)-
      x*exp(-beta*x))/((1/(exp(beta)-1))+exp(-beta*x)) ))^2)*
    ((1/beta)/( (1/(exp(beta)-1)) + exp(-beta*x)))
  mean(y)
}
### FD type two
fishinf.FD2_2=function(beta,n){
  x=rexp(n,rate=beta)
  y=(((1/beta)-x+(x*exp(-beta*x))/(exp(-beta*x)+1))^2)*
    (1)/((exp(-beta*x)+1)*log(2,base=exp(1)))
  mean(y)
}

#example
fishinf.FD1_2(beta = 0.25, n = 1e6)
fishinf.FD2_2(beta = 0.25, n = 1e6)

##### entropy R integral ######
### FD type one
ent.FD1_1=function(beta){
  f=function(x,beta){
  -((exp(-beta*x))/((1)/(exp(beta)-1)+exp(-beta*x)))*
      (-beta*x-log((1/(exp(beta)-1))+exp(-beta*x),base=exp(1)))
  }
  integrate(f,beta,lower=0,upper=Inf)
}
### FD type two
ent.FD2_1=function(beta){
  f=function(x,beta){
    -((beta*exp(-beta*x))/(log(2,base=exp(1))*(1+exp(-beta*x))))*
      (log(beta,base=exp(1))-beta*x-log(log(2,base=exp(1)),exp(1))-
         log(exp(-beta*x)+1,base=exp(1)))
  }
  integrate(f,beta,lower=0,upper=Inf)
}

#example
ent.FD1_1(beta = 0.25)
ent.FD2_1(beta = 0.25)

##### entropy MCMC #####
### FD type one
ent.FD1_2=function(beta,n){
  x=rexp(n,rate=beta)
  y=((1/beta)/((1)/(exp(beta)-1)+exp(-beta*x)))*
    (beta*x+log( (1/(exp(beta)-1))+exp(-beta*x),base=exp(1)))
  mean(y)
}
### FD type two
ent.FD2_2=function(beta,n){
  x=rexp(n,rate=beta)
  y=((1)/(log(2,base=exp(1))*(1+exp(-beta*x))))*
    (-log(beta,base=exp(1))+beta*x+log(log(2,base=exp(1)),exp(1))+
       log(exp(-beta*x)+1,base=exp(1)))
  mean(y)
}

#example
ent.FD1_2(beta = 0.25, n = 1e6)
ent.FD2_2(beta = 0.25, n = 1e6)

##### generate random sample #####
### FD type one
rs.FD1=function(beta,a,n){
  u=runif(n)
  -(log(exp(-beta*u)*((exp(-beta*a))/(exp(beta)-1)+
      exp(-beta*a))-(exp(-beta*a))/(exp(beta)-1),base=exp(1)))/
    (beta)
}
### FD type two
rs.FD2=function(beta,a,n){
  u=runif(n)
  -(log(exp(-u*log(1+exp(-beta*a),base=exp(1)))*
      (1+exp(-beta*a))-1,base=exp(1)))/(beta)
}

# n is the size of the sample

#example
rs.FD1(beta = 0.25, a = 0, n = 5)
rs.FD2(beta = 0.25, a = 0, n = 5)

##### CDF ######
### FD type one
cdf.FD1=function(x,beta,a=0){
  (-1/beta)*log(((exp(-beta*a)/(exp(beta)-1))+exp(-beta*x))/
        ((exp(-beta*a)/(exp(beta)-1))+exp(-beta*a)),base=exp(1))
}
### FD type two
cdf.FD2=function(x,beta,a=0){
  -(1/log(1+exp(-beta*a),base=exp(1)))*
    log( (1+exp(-beta*x))/(1+exp(-beta*a)),base=exp(1))
}

#example
cdf.FD1(x = 1, beta = 0.25, a = 0)
cdf.FD2(x = 1, beta = 0.25, a = 0)

##### KL divergence #####
### FD1 vs FD2   D(f_FD_1 || f_FD_2)
kl.FD1FD2=function(beta1, beta2, n){
  
  x1=rexp(n,beta1)
  y1=((1/beta1)/((1/(exp(beta1)-1))+exp(-beta1*x1)))*
    log(( exp(-beta1*x1)*(exp(-beta2*x1)+1)*log(2,base=exp(1)))/
      (beta2*exp(-beta2*x1)*((1/(exp(beta1)-1))+exp(-beta1*x1)))
      ,base=exp(1))
  fd1fd2=mean(y1)
  
  x2=rexp(n,beta2)
  y2=(1/((exp(-beta2*x2)+1)*log(2,base=exp(1))))*
    log((beta2*exp(-beta2*x2)*((1/(exp(beta1)-1))+exp(-beta1*x2)))/
      (exp(-beta1*x2)*(exp(-beta2*x2)+1)*log(2,base=exp(1)))
      ,base=exp(1))
  fd2fd1=mean(y2)
  
  dbar=(fd1fd2+fd2fd1)/2
  
  list(FD1FD2=fd1fd2 , FD2FD1=fd2fd1 , Dbar=dbar)
}

# beta1 is the parameter of FD1 and beta2 is the parameter
# of FD2. n is the size of the sample.

#example
kl.FD1FD2(beta1 = 0.25, beta2 = 1, n = 1e6)

##### log likelihood #####
### FD type one
lnlike.FD1=function(samp,beta){
  -beta*sum(samp)-sum(log((1/(exp(beta)-1))+
              exp(-beta*samp),base=exp(1)))
}
### FD type two
lnlike.FD2=function(samp,beta){
  length(samp)*log(beta,base=exp(1))-beta*sum(samp)-
    length(samp)*log(log(2,base=exp(1)),base=exp(1))-
    sum(log(exp(-beta*samp)+1,base=exp(1)))
}

# samp is the vector of sample

#example
lnlike.FD1(samp = rs.FD1(beta = 0.25, a = 0, n = 5), beta = 0.25)
lnlike.FD2(samp = rs.FD2(beta = 0.25, a = 0, n = 5), beta = 0.25)

### median ####
### FD type one
med.FD1=function(beta,a){
  -(1/beta)*log(exp(-beta/2)*((exp(-beta*a)/
      (exp(beta)-1))+exp(-beta*a))-(exp(-beta*a)/
                  (exp(beta)-1)),base=exp(1))
}
### FD type two
med.FD2=function(beta,a){
  -(1/beta)*log(exp(-(1/2)*log(1+exp(-beta*a),base=exp(1)))*
                  (exp(-beta*a)+1)-1,base=exp(1))
}

# example
med.FD1(beta = 0.25, a = 0)
med.FD2(beta = 0.25, a = 0)

##### GLRT #####
### FD type one
glrt.FD1=function(samp,beta0){
  mle=optimize(lnlike.FD1,samp=samp,interval = c(0,50),
               maximum = T)$maximum
  test.stat=2*(lnlike.FD1(samp,mle) - lnlike.FD1(samp,beta0))
  reject.area=qchisq(0.95,1)
  p.value=1-pchisq(test.stat,1)
  res=if(test.stat>reject.area) 'reject' else 'accepted'
  list(test.statistic = test.stat, rejection_area = reject.area,
       p_value=p.value, result=res)
}
### FD type two
glrt.FD2=function(samp,beta0){
  mle=optimize(lnlike.FD2,samp=samp,interval = c(0,50),
               maximum = T)$maximum
  test.stat=2*(lnlike.FD2(samp,mle)-lnlike.FD2(samp,beta0))
  reject.area=qchisq(0.95,1)
  p.value=1-pchisq(test.stat,1)
  res=if(test.stat>reject.area) 'reject' else 'accepted'
  list(test.statistic = test.stat, rejection_area = reject.area,
       p_value=p.value, result=res)
}

# samp is the vector of sample, beta0 is the value in null
# hypothesis.

# example
glrt.FD1(samp = rs.FD1(beta = 0.25, a = 0, n = 100), beta0 = 0.25)
glrt.FD2(samp = rs.FD2(beta = 0.25, a = 0, n = 100), beta0 = 0.25)

##### Wald test #####
### FD type one
waldt.FD1=function(samp,beta0){
  mle=optimize(lnlike.FD1,samp=samp,interval = c(0,50),
               maximum = T)$maximum
  var.mle=length(samp)*fishinf.FD1_1(mle)$value
  CI=c(lower=mle-qnorm(0.975)*(1/sqrt(var.mle)),
       upper=mle+qnorm(0.975)*(1/sqrt(var.mle)))
  test.stat=(mle-beta0)/(1/sqrt(var.mle))
  res= if(abs(test.stat)>qnorm(0.975)) 'reject' else 'accept'
  p.value=2*(1-pnorm(abs(test.stat)))
  
  list('maximum likelihood' = mle , 'variance mle' = 1/var.mle,
       'confidence interval' = CI , 'test statistic' = test.stat,
       'result' = res, 'p_value' = p.value)
}

### FD type two
waldt.FD2=function(samp,beta0){
  mle=optimize(lnlike.FD2,samp=samp,interval = c(0,50),
               maximum = T)$maximum
  var.mle=length(samp)*fishinf.FD2_1(mle)$value
  CI=c(lower=mle-qnorm(0.975)*(1/sqrt(var.mle)),
       upper=mle+qnorm(0.975)*(1/sqrt(var.mle)))
  test.stat=(mle-beta0)/(1/sqrt(var.mle))
  res= if(abs(test.stat)>qnorm(0.975)) 'reject' else 'accept'
  p.value=2*(1-pnorm(abs(test.stat)))
  list('maximum likelihood' = mle , 'variance mle' = 1/var.mle,
       'confidence interval' = CI , 'test statistic' = test.stat,
       'result' = res, 'p_value' = p.value)
}

# example
waldt.FD1(samp = rs.FD1(beta = 0.25, a = 0, n = 100), beta0 = 0.25)
waldt.FD2(samp = rs.FD2(beta = 0.25, a = 0, n = 100), beta0 = 0.25)

### log likelihood prime ####
### FD type one
lnlikeprim.FD1=function(samp,beta){
  -sum(samp)-sum( ( (-exp(beta)/(exp(beta)-1)^2)-
      samp*exp(-beta*samp) )/( (1/(exp(beta)-1))+
                          exp(-beta*samp) ) )
}
### FD type two
lnlikeprim.FD2=function(samp,beta){
  length(samp)/beta - sum(samp) + sum((samp*exp(-beta*samp))/
                    (exp(-beta*samp)+1))
}

# example
lnlikeprim.FD1(samp = rs.FD1(beta = 0.25, a = 0, n = 10), beta = 0.25)
lnlikeprim.FD2(samp = rs.FD2(beta = 0.25, a = 0, n = 10), beta = 0.25)

### score test ####
### FD type one
score.FD1=function(samp,beta0){
  s=(lnlikeprim.FD1(samp,beta = beta0)^2)/
    (length(samp)*fishinf.FD1_1(beta=beta0)$value)
  rej=qchisq(0.95,1)
  res= if(s>rej) 'reject' else 'accept'
  p.value=1-pchisq(s, 1)
  
  list('test statistic' = s, 'rejection area' = rej,
       'result' = res, 'p_value' = p.value)
}
### FD type two
score.FD2=function(samp,beta0){
  s=(lnlikeprim.FD2(samp,beta = beta0)^2)/
    (length(samp) * fishinf.FD2_1(beta=beta0)$value)
  rej=qchisq(0.95,1)
  res= if(s>rej) 'reject' else 'accept'
  p.value=1-pchisq(s, 1)
  
  list('test statistic' = s, 'rejection area' = rej,
       'result' = res, 'p_value' = p.value)
}

# example
score.FD1(samp = rs.FD1(beta = 0.25, a = 0, n = 100), beta0 = 0.25)
score.FD2(samp = rs.FD2(beta = 0.25, a = 0, n = 100), beta0 = 0.25)


##### thank you for your attention.