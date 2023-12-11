


genModWish2=function(dat, M = 2000, b0 = .03){
  if (mean(dat$cond %in% 1:2)<1) stop("Conditions must be 1 and 2")
  dat$sub=as.integer(as.factor(dat$sub))
  I=max(dat$sub)
  J=max(dat$task)
  N=dim(dat)[1]
  K=table(dat$sub,dat$task,dat$cond)
  mn=tapply(dat$rt,list(dat$sub,dat$task,dat$cond),mean)
  sd=tapply(dat$rt,list(dat$sub,dat$task,dat$cond),sd)
  x=tapply(dat$cond,list(dat$sub,dat$task,dat$cond),mean)-1
  
  theta=alpha=array(dim=c(M,I,J),0)
  s2=1:M
  muTheta=matrix(nrow=M,ncol=J)
  s2[1]=.25^2
  s2Alpha=1
  muAlpha=.8
  muTheta[1,]=rep(0,J)
  
  meanMuTheta=rep(.05,J)
  precMuTheta=diag(J)/.1^2
  B=array(dim=c(M,J,J))
  B[1,,]=diag(J)/.025^2
  XtX=K[,,2]
  
  a0=J
  for (m in 2:M){
    for (i in 1:I) {
      #alpha
      c=apply(K[i,,]*(mn[i,,]-x[i,,]*theta[m-1,i,]),1,sum)/s2[m-1]+muAlpha/s2Alpha
      v=1/(apply(K[i,,],1,sum)/s2[m-1]+1/s2Alpha)
      alpha[m,i,]=rnorm(J,c*v,sqrt(v))
      #theta
      Xty=K[i,,2]*(mn[i,,2]-alpha[m,i,])
      c=Xty/s2[m-1]+B[m-1,,]%*%muTheta[m-1,]
      v=solve(diag(XtX[i,])/s2[m-1]+B[m-1,,])
      theta[m,i,]=rmvnorm(1,v%*%c,v)
    }
    #s2
    scale=sum((K-1)*sd^2+K*((mn-outer(alpha[m,,],c(1,1))-x*outer(theta[m,,],c(1,1)))
                            ^2))/2+.5
    s2[m]=rinvgamma(1,shape=(N+1)/2,scale=scale)
    #muTheta
    v=solve(I*B[m-1,,]+precMuTheta)
    c=I*B[m-1,,]%*%apply(theta[m,,],2,mean)+precMuTheta%*%meanMuTheta
    muTheta[m,]=rmvnorm(1,v%*%c,v)
    #B
    err=t(t(theta[m,,])-muTheta[m,])
    SSE=crossprod(err)
    scale=solve(diag(rep(b0^2,J))  + SSE)
    B[m,,]=rWishart(1,a0+I,scale)
  }
  return(list(alpha=alpha,theta=theta,B=B))
}


modIWJags = "
model{
  # Priors
  for (i in 1:I){
    alpha[i, 1:J] ~ dmnorm(mu, pPsi)
    theta[i, 1:J] ~ dmnorm(nu, pDel)
  }
  
  for (j in 1:J){
    mu[j] ~ dnorm(.8, pow(.3, -2))
    nu[j] ~ dnorm(0, pow(.1, -2))
  }
  
  pPsi ~ dwish(diagJ*tuneA^2, J+1)
  pDel ~ dwish(diagJ*tuneT^2, J+1)
  pSig ~ dgamma(.5,.5)
  
  
  # Likelihood
  for (n in 1:N){
    center[n] = alpha[sub[n], task[n]] + (cond[n] - 1.5) * theta[sub[n], task[n]]
    y[n] ~ dnorm(center[n], pSig)
  }
}
"