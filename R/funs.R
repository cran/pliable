#' Fit a  pliable lasso model over a path of regularization values
#' @param x  n by p matrix of predictors
#' @param z n by nz matrix of modifying variables. The elements of z  may represent quantitative or categorical variables, or a mixture of the two.
#' Categorical varables should be coded by 0-1 dummy variables: for a k-level variable, one can use either k or k-1  dummy variables.
#' @param y n-vector of responses. The x and z variables are centered in the function. We recommmend that x and z also be standardized before the call
#'   via x<-scale(x,T,T)/sqrt(nrow(x)-1);   z<-scale(z,T,T)/sqrt(nrow(z)-1)
#' @param lambda decreasing sequence of lam values desired. Default NULL.  If NULL, computed in function
#' @param nlambda number of lambda values desired (default 50).
#' @param family response type- either "gaussian", "binomial". In the binomial case, y should be 0s and 1s.
#'    family "cox" (Prop Hazards model for right-censored data) will be implemented if there is popular demand
#' @param lambda.min.ratio the smallest value for lambda , as a fraction of
#'     lambda.max, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero).
#'      Default is 0.001 if n>p, and 0.01 if n< p.
#'      Along with "screen" (see below),  this  parameter is a main parameter
#'      controlling the runtime of the algorithm, With a  smaller lambda.min.ratio,  more interactions will typically be entered.
#'       If the algorithm stalls or stops near the end of the path, trying increasing lambda.min.ratio
#'        On the other hand, if a more complex fit is desired, or  the cross-validation curve from cv.pliable is
#'      still going down at the end of the path,
#'        consider decreasing lambda.min.ratio
#' @param alpha mixing parameter- default 0.5
#' @param w n-vector of observation wts (Default: all ones). Not currently used with Cox family.
#' @param thr convergence threshold  for estimation of beta and joint estimation of (beta,theta).
#'      Default 1e-5
#' @param  penalty.factor vector of penalty factors, one for each X  predictors. Default: a vector of ones
#' @param tt  initial factor for Generalized grad algorithm for joint estimation of main effects and
#'    interaction parameters- default 0.1/mean(x^2)
#' @param maxit maximum number of iterations in loop for one lambda. Default 10000
#' @param mxthit maximum number of theta solution iterations. Default 1000
#' @param mxkbt maximum number of backtracking iterations. Default 100
#' @param mxth maximum internal theta storage. Default 1000
#' @param kpmax maximum dimension of main effect storage. Default 100000
#' @param kthmax maximum dimension of interaction storage. Default 100000
#' @param verbose  Should information  be printed along the way?  Default FALSE
#' @param maxinter Maximum number of interactions allowed in the model. Default 500
#' @param zlinear Should an (unregularized) linear term for z be included in the model? Default TRUE.
#' @param screen Screening quantile for features: for example a value screen=0.7  means that we keep only those features with univariate t-stats above the 0.7 quantile of all features.
#'  Default is NULL, meaning that all features are kept
#'
#'
#' @details
#' This function fits a pliable lasso model over a sequence of lambda (regularization) values.
#' The inputs are a feature matrix x, response vector y
#' and a matrix of modifying variables z.  The pliable lasso is a generalization of the lasso in which the weights multiplying the features x
#' are allowed to vary as a function of a set of modifying variables z. The model has the form:
#'
#'  yhat=beta0+ z\%*\% betaz+ rowSums( x[,j]*beta[j] + x[,j]*z\%*\% theta[j,])
#'
#' The parameters are beta0 (scalar), betaz (nz-vector), beta (p-vector) and theta (p by nz).

#'  The  (convex) optimization problem is
#'
#'    minimize_(beta0,betaz,beta,theta)  (1/(2n))* sum_i ((y_i-yhat_i)^2) +(1-alpha)*lambda* sum_j (||(beta_j,theta_j)||_2 +||theta_j)|| +alpha*lambda* sum_j||theta_j||_1
#'
#'  A blockwise coordinate descent procedure is employed
#' for the optimization.
#'
#' @return  A trained  pliable object with the following components:
#'
#' param:  values of lambda used
#'
#'  beta: matrix of estimated beta coefficients  ncol(x) by length(lambda).
#'
#'  theta: array of estimated theta coefficients  ncol(x) by ncol(z) by length(lambda).
#'
#' betaz: estimated (main effect)  coefficients for z  ncol(z) by length(lambda).
#'
#' xbar:  rowMeans of x
#'
#' zbar:  rowMeans of z
#'
#' w:     observation weights used
#'
#'  args: list of arguments in the original call
#'
#' nulldev:   null deviance for model
#'
#' dev.ratio:  deviance/null.dev for each model in path
#'
#' nbeta:    number of nonzero betas for each model in path
#'
#' nbeta.with.int:  number of  betas with some nonzero theta values, for each model in path
#'
#' ntheta: number of nonzero theta for each model in path
#'
#' df: rough degrees of freedom for each model in path  (=nbeta)
#'
#' istor: for internal use
#'
#' fstor:  for internal use
#'
#'  thstor: for internal use
#'
#' jerr: for internal use
#'
#' call:  the call made to pliable
#'
#'
#'

#'
#' @examples
#' # Train a pliable lasso model- Gaussian case
#' #  Generate some data
#' set.seed(9334)
#' n = 20; p = 3 ;nz=3
#' x = matrix(rnorm(n*p), n, p)
#' mx=colMeans(x)
#' sx=sqrt(apply(x,2,var))
#' x=scale(x,mx,sx)
#' z =matrix(rnorm(n*nz),n,nz)
#' mz=colMeans(z)
#' sz=sqrt(apply(z,2,var))
#' z=scale(z,mz,sz)
#' y =4*x[,1] +5*x[,1]*z[,3]+ 3*rnorm(n)
#'  # Fit model
#'  fit = pliable(x,z,y)
#'  fit   #examine results
#' plot(fit)  #plot path
#'
#' #Estimate main tuning parameter lambda by cross-validation
#'  #NOT RUN
#'  # cvfit=cv.pliable(fit,x,z,y,nfolds=5) # returns lambda.min and lambda.1se,
#'     #  the minimizing and 1se rules
#'  # plot(cvfit)
#'
#' # Predict using the fitted model
#' #ntest=500
#' #xtest = matrix(rnorm(ntest*p),ntest,p)
#' #xtest=scale(xtest,mx,sx)
#' #ztest =matrix(rnorm(ntest*nz),ntest,nz)
#' #ztest=scale(ztest,mz,sz)
#' #ytest = 4*xtest[,1] +5*xtest[,1]*ztest[,3]+ 3*rnorm(ntest)
#' #pred= predict(fit,xtest,ztest,lambda=cvfit$lambda.min)
#' #plot(ytest,pred)
#'
#' #Binary outcome example
#' n = 20 ; p = 3 ;nz=3
#' x = matrix(rnorm(n*p), n, p)
#' mx=colMeans(x)
#' sx=sqrt(apply(x,2,var))
#' x=scale(x,mx,sx)
#' z =matrix(rnorm(n*nz),n,nz)
#' mz=colMeans(z)
#' sz=sqrt(apply(z,2,var))
#' z=scale(z,mz,sz)
#' y =4*x[,1] +5*x[,1]*z[,3]+ 3*rnorm(n)
#' y=1*(y>0)
#'  fit = pliable(x,z,y,family="binomial")
#'
#'
#'
#' # Example where z is not observed in the test set, but predicted
#'  #        from a supervised learning algorithm
#'
#'
#' n = 20; p = 3 ;nz=3
#' x = matrix(rnorm(n*p), n, p)
#' mx=colMeans(x)
#' sx=sqrt(apply(x,2,var))
#' x=scale(x,mx,sx)
#' z =matrix(rnorm(n*nz),n,nz)
#' mz=colMeans(z)
#' sz=sqrt(apply(z,2,var))
#' z=scale(z,mz,sz)
#' y =4*x[,1] +5*x[,1]*z[,3]+ 3*rnorm(n)
#'
#' fit = pliable(x,z,y)
#' # predict z  from x; here we use glmnet, but any other supervised method can be used
#'
#'
#' zfit=cv.glmnet(x,z,family="mgaussian")
#'
#' # Predict using the fitted model
# # NOT RUN
#' # ntest=100
#' #xtest =matrix(rnorm(ntest*nz),ntest,p)
#' #xtest=scale(xtest,mx,sx)
#' #ztest =predict(zfit,xtest,s=cvfit$lambda.min)[,,1]
#' #ytest = 4*xtest[,1] +5*xtest[,1]*ztest[,3]+ 3*rnorm(ntest)
#'
#' #pred= predict(fit,xtest,ztest)
#' #plot(ytest,pred[,27])  # Just used 27th path member, for illustration; see cv.pliable for
#' #  details of how to cross-validate in this setting
#'
#'
#'
#' @seealso
#' cv.pliable, predict.pliable, plot.pliable, plot.cv.pliable
#'
#' @export
#'
#'

pliable <- function(x,z,y,lambda=NULL,nlambda=50,alpha=.5,family=c("gaussian","binomial"),
               lambda.min.ratio=NULL,
               w=NULL,thr=1e-5,tt=NULL,penalty.factor=rep(1,ncol(x)), maxit=10000,mxthit=100,mxkbt=100,
                mxth=1000,
                kpmax=100000,kthmax=100000,verbose=FALSE, maxinter=500, zlinear=TRUE, screen=NULL){
    this.call = match.call()

 #   debug=F
 #   if(debug){
 #       lambda=NULL
 #   nlambda=50
 #   alpha=.5
  #  lambda.min.ratio=NULL
 #   w=rep(1,length(y))
  #  thr=1e-5
  #  tt=NULL
  #  maxit=10000
 #  mxthit=100
 #   mxkbt=100
 #   mxth=1000
 #   kpmax=100000
 #   kthmax=100000
 #  verbose=TRUE
 #   maxinter=500
 #   zlinear=TRUE
 #   screen=NULL
 #   }
    SMALL=1e-8
    family=match.arg(family)
    if(family=="cox" & !is.null(w)) cat("wts ignored for Cox family",fill=T)
        if(is.null(w)) w=rep(1,length(y))
     if(family=="gaussian") errfun=errfun.gaussian
    if(family=="binomial") errfun=errfun.binomial
     if(family=="cox")  errfun=function(yhat,coxinfo) {devc(no,coxinfo$kq, coxinfo$ddq, yhat, coxinfo$iriskq,status)}

    if(!is.matrix(x)) x=matrix(x,ncol=1)
    if(!is.matrix(z)) z=matrix(z,ncol=1)

    if(family=="cox"){
        # order data by y
        o=order(y-.0001*status)
        x= x[o,]
        z=z[o,]
        y=y[o]
        status=status[o]
        w=w[o]
    }
           #   browser()
    coxinfo=NULL
    if(family=="cox") coxinfo=initcx(x,y,status)

    vz=apply(z,2,var)
    if(any(vz< SMALL)){
          mess=which(vz<SMALL)[1];for(i in which(vz<SMALL)[-1] ){
                   mess=paste(mess,",",i,sep="")
                                  }
 #        browser()
                   stop(c("Columns ", mess, " of z are constant"))
    }

    xbar=colMeans(x)
    zbar=colMeans(z)

    x=scale(x,xbar,F)
    z=scale(z,zbar,F)

    orig.x=x
    if(!is.null(screen)){
         ni.original=ncol(x)
         if(family=="gaussian") ttv=quantitative.func(t(x),y)$tt
        if(family=="binomial") ttv=ttest.func(t(x),y)$tt
        kee=which(abs(ttv)>= quantile(abs(ttv),screen))
        x=x[,kee]
    }

  #  browser()

    no=nrow(x)
ni=ncol(x)
nz=ncol(z)
    nt=ni+nz
    y=as.vector(y)




 if(is.null(lambda.min.ratio)){
     if( (family=="gaussian" | family=="cox") & no>ni) lambda.min.ratio=.001
        if( (family=="gaussian" | family=="cox") & no<ni) lambda.min.ratio=.01
        if(family=="binomial" & no>ni) lambda.min.ratio=.001
        if(family=="binomial" & no<ni) lambda.min.ratio=.01
    }

## error checking of input args

 if (alpha > 1) {
        warning("alpha >1; set to 1")
        alpha = 1
 }

    if(any(penalty.factor<0)) stop("Penalty factors should be ge 0")
  if(any(penalty.factor>0)) penalty.factor=ni*penalty.factor/sum(penalty.factor)

    if (alpha < 0) {
        warning("alpha<0; set to 0")
        alpha = 0
    }
    if(!is.null(lambda)) nlambda=length(lambda)

    if(!is.null(lambda)) {if(any(lambda<0)) stop("lambda must be non-negative")}


    if(family=="binomial" & (sum(y==0)+sum(y==1))!=no) stop("y must be 0s and 1s for binomial family")
    if(family=="cox"){
        if( is.null(status)) stop("with family='cox'  need to provide status argument")
         if( any(y<=0)) stop("outcome y must be non-negative with  family='cox'")
        if((sum(status==1)+sum(status==0))!=length(y))  stop("status variable must be 1 or 0")
        }
###

    if(is.null(lambda)){
        r=y
        if(zlinear) {
            r=lsfit(z,y)$res
            if(sum(r^2)< 0.0001) cat("linear model in z overfit; consider rerunning with zlinear=F",fill=T)
            }
        if(family=="gaussian" | family=="binomial") lammax=max(abs(t(x)%*%r)/length(r))/(1-alpha)

        if(family=="cox") {

                          junk2=calcz(x,coxinfo$kq,rep(0,no),status,coxinfo$iriskq,coxinfo$ddq,rep(0,ncol(x)))
            ll <- as.vector(junk2$wz)
              l <- as.vector(junk2$grad)
             r = - l/ll
             if(zlinear){
               wt=ll;wt[wt<0]=0
          zcoef=lsfit(z,r,intercept=F,wt=wt)$coef
          }

        lammax=max(abs(t(x) %*% (-ll * r))/no)/(1-alpha)
      }
        lambda= exp(seq(log(lammax),log(lammax*lambda.min.ratio),length=nlambda))

        }



    yhat=rep(0,no)

    ybar=0

    if(family=="gaussian"){
    ybar=mean(y)
    y=y-ybar
    }


nlambda=length(lambda)
xx=cbind(x,z)

    ulam=lambda

    if(is.null(tt)) tt=.1/mean(x^2)
    if(verbose) {cat(c("initial value of backtrack parameter tt=",tt),fill=T);cat("",fill=T)}

    iverbose=1*verbose
    kbos=iverbose
    linenter=1*zlinear
    mlam=0
mode(xx)="double"
mode(y)="double"
mode(w)="double"
mode(alpha)="double"
mode(nlambda)="integer"
    mode(ulam)="double"
    mode(kbos)="integer"
    mode(maxinter)="integer"

    mode(linenter)="integer"
    mode(thr)="double"
    mode(penalty.factor)="double"
    mode(tt)="double"
    mode(maxit)="integer"
    mode(mxthit)="integer"
     mode(mxkbt)="integer"
mode(mxth)="integer"
mode(kpmax)="integer"
mode(kthmax)="integer"

mode(mlam)="integer"

 #   if(!is.loaded("pliable"))  dyn.load("/Users/tibs/dropbox/PAPERS/weightree2/jerry/pliable2.so")

  #  if(family=="binomial" & !is.loaded("logpliable"))  dyn.load("/Users/tibs/dropbox/PAPERS/weightree2/jerry/logpliable.so")



    if(family=="gaussian"){
out=.Fortran("pliable",
             no,ni,nz,xx,y,w,alpha,nlambda,ulam,kbos,maxinter,linenter,thr,penalty.factor,maxit,tt,mxthit, mxkbt,mxth,kpmax,kthmax,mlam,
   a0=double(nlambda),
   lp=integer(2*nlambda),
   istor=integer(2*kpmax),
   fstor=double(kpmax),
   thstor=double(kpmax),
   jerr=integer(1),
   PACKAGE="pliable"
)
    }

    if(family=="binomial"){
out=.Fortran("logpliable",
            no,ni,nz,xx,y,w,alpha,nlambda,ulam,kbos,maxinter,linenter,thr,penalty.factor,maxit,tt,mxthit, mxkbt,mxth,kpmax,kthmax,mlam,
       a0=double(nlambda),
   lp=integer(2*nlambda),
   istor=integer(2*kpmax),
   fstor=double(kpmax),
   thstor=double(kpmax),
   jerr=integer(1),
     PACKAGE="pliable"
)
    }



    if(family=="cox"){
        mode(status)="integer"
out=.Fortran("coxpliable",
             no,ni,nz,xx,y,status,w,alpha,nlambda,ulam,kbos,maxinter,linenter,
           thr,maxit,tt,mxthit, mxkbt,mxth,kpmax,kthmax,mlam,
       a0=double(nlambda),
   lp=integer(2*nlambda),
   istor=integer(2*kpmax),
   fstor=double(kpmax),
   thstor=double(kpmax),
   jerr=integer(1),
   PACKAGE="pliable"
)
    }

    out[1:17]=NULL  #remove blank stuff for return

#
    if(out$jerr<0){
        cat(c("jerr=",out$jerr),fill=T)
       cat(c("max iters exceeded, jerr=",out$jerr),fill=T)
    }
     if(out$jerr>0){
        cat(c("jerr=",out$jerr),fill=T)
        stop(c("Storage error in pliable, jerr=",out$jerr))
        }

     out$istor=matrix(out$istor,nrow=2)
    out$lp=matrix(out$lp,nrow=2)
    out$lambda=lambda
     out$args=list(nlambda=nlambda,alpha=alpha,lambda.min.ratio=lambda.min.ratio,w=w,thr=thr,penalty.factor=penalty.factor,tt=tt,maxit=maxit,mxthit=mxthit,mxkbt=mxkbt,mxth=mxth,kpmax=kpmax,kthmax=kthmax,verbose=verbose,screen=screen)

    beta=matrix(0,nrow=ncol(x),ncol=nlambda)
    betaz=matrix(0,nrow=ncol(z),ncol=nlambda)
    theta=array(0,c(ncol(x),ncol(z),nlambda))
    nz.actual=nz*(zlinear)

    for(k in 1:nlambda){
  #      cat(k)
 #       if(k==18) browser()
        if(ncol(out$lp)>0 ){
            if(out$lp[1,k]>0){

           junk=modsoln(out,ni,nz,k)
 #browser()
        nb=length(junk$iv)-nz.actual
        kv=junk$kv
        if(nb>0)  beta[junk$iv[1:nb],k]=junk$av[1:nb]
        zlist=(nb+1):length(junk$iv)
        if(zlinear){
        iz=junk$iv[zlist]-ncol(x)
        betaz[iz,k]=junk$av[zlist]
        }
#   kv = number of non-zero coefficients
#   iv(kv) = identies of non-zero coefficients
#   av(kv) = corresponding non-zero values
#   it(kv) = interaction pointer for each non-zero coefficient:
#      it(j)=0 => no interaction for iv(j) coefficient
#      it(j)>0 => interaction coefficients (thetas) in th(1:nz,it(j))
#   kz = number of interacting coefficient sets for this solution
#   th(nz,kz) = all interaction coefficient sets for this solution

        th=t(matrix(junk$th,nrow=nz))
      #      iv=rev(rev(junk$iv)[-(1:nz)]) #remove last nz components corr to z vars
      #      it=rev(rev(junk$it)[-(1:nz)])
            inds=junk$iv[junk$it>0]

         #   inds=inds[inds<=ni]
        if(length(inds)>0)   theta[inds, ,k]=th[1:length(inds),]
    }}}

   out$family=family
    out$beta=beta
    out$df=colSums(beta!=0)
    out$theta=theta
   out$betaz=betaz
    out$xbar=xbar
    out$zbar=zbar
    out$ybar=ybar

    out$w=w



    if(!is.null(screen)){
        tbeta=matrix(0,ni.original,ncol(out$beta))
        tbeta[kee,]=out$beta
        out$beta=tbeta
         ttheta=array(0,c(ni.original,dim(out$theta)[2],dim(out$theta)[3]))
        ttheta[kee,,]=out$theta
        out$theta=ttheta
        }

    # have to uncenter x and z below, to make it work with predict function
    if(family=="gaussian") {

          pred=predict.pliable(out,scale(orig.x,-xbar,FALSE), scale(z,-zbar,FALSE),type="response")
        dev=colMeans( errfun(y+ybar,pred,w))
        out$nulldev=mean(errfun(y+ybar,ybar,w))
    }
     if( family=="binomial") {

          pred=predict.pliable(out,scale(orig.x,-xbar,FALSE), scale(z,-zbar,FALSE),type="response")
        dev=colMeans( errfun(y,pred,w))
         out$nulldev=mean(errfun(y,mean(y),w))
        }
    if(family=="cox"){
          pred=predict.pliable(out,scale(orig.x,-xbar,FALSE), scale(z,-zbar,FALSE),type="link")
     dev=apply(pred,2,errfun,coxinfo)
     out$nulldev=errfun(rep(0,length(y)),coxinfo)
     }
  #  dev=dev[!is.na(dev)]
  #  out$nulldev=mean( (y-mean(y))^2)

    out$dev.ratio= dev/out$nulldev
    out$dev.ratio=out$dev.ratio[!is.na(out$dev.ratio)]
    oo=rev(colSums(abs(out$beta!=0)))  # find zeros at end of path and remove them
    o2=which(cumsum(oo)==0)
    if(length(o2)>0){
     ooo=max(o2)
     ooo=max(ooo,0)
     nonzero=1:(length(out$dev.ratio)-ooo)
     out$dev.ratio=out$dev.ratio[nonzero]
     out$beta=out$beta[,nonzero]
     out$theta=out$theta[,,nonzero]
     out$betaz=out$betaz[,nonzero]
     out$lambda=out$lambda[nonzero]
    }
     if(any(diff(out$dev.ratio)>.01)) cat(c("Caution: deviance increasing along path; rerun with verbose=TRUE to
see the value of backtrack parameter tt, and then try rerunning pliable with tt= half its current value of ",tt))
    out$nbeta=colSums(out$beta!=0)
    out$nbeta.with.int=apply(apply(out$theta!=0,c(1,3),sum),2,sum)
    out$ntheta=apply(out$theta!=0,c(3),sum)
    out$coxinfo=coxinfo
    out$pred=pred
out$call=this.call
    class(out)="pliable"
return(out)
}



tn=function(x) sqrt(sum(x*x))


modsoln=function(out,ni,nz,k){

# browser()
istor=out$istor
fstor=out$fstor
thstor=out$thstor

lpp=out$lp[,k]
mode(lpp)="integer"
mode(istor)="integer"
mode(fstor)="double"
mode(nz)="integer"
mode(thstor)="double"

mxth=out$args$mxth



out2=.Fortran("modsoln",
              nz,lpp,istor,fstor,thstor,kv=integer(1) ,iv=integer(500000),av=double(500000),it=integer(500000),kz=integer(1),
              th=double(nz*mxth))
       #make sure dims of iv,av,it ar big enough!

#if(k==25) browser()
iv=out2$iv[1:out2$kv]
av=out2$av[1:out2$kv]
it=out2$it[1:out2$kv]
th=out2$th


out=list(kv=out2$kv,iv=iv,av=av,it=it,kz=out2$kz,th=th)

return(out)
}

#' Compute predicted values from a fitted pliable  object

#' Make predictions from a fitted pliable lasso model
#' @param object object returned from a call to pliable
#' @param x  n by p matrix of predictors
#' @param z  n by nz matrix of modifying variables. These may be observed or the predictions from a supervised learning
#'   algorithm that predicts z from test features x  and possibly other features. See example below

#'
#' @param type Returns either the fitted values with type="link" or  "response"  or the parameter estimates when type="coefficients" . Type "link" gives the linear
#'          predictors for binomial,  or
#'          cox models; for gaussian models it gives the fitted
#'          values. Type "response" gives the fitted probabilities for
#'          binomial,
#'          and the fitted relative-risk for cox; for gaussian
#'          type "response" is equivalent to type "link".
#' @param lambda  values of lambda at whcih predictions are desired. If NULL (default), the path of lambda values from the fitted model.
#'   are used. If lambda is  not NULL,  the predictions are made at the closest values to lambda in the lambda path  from  the fitted model;
#' @param verbose  Should information should be printed along the way?  Default FALSE.
#' @param ... Further arguments (not used)
#'  @return  predicted values
#'
#' @examples
#' # Train a pliable lasso model
#' n = 20; p = 3 ;nz=3
#' x = matrix(rnorm(n*p), n, p)
#' z =matrix(rnorm(n*nz),n,nz)
#' y = x[,1] +x[,1]*z[,3]+ rnorm(n)
#'   fit = pliable(x,z,y)
#'  # plot coefficient profiles, indicating z-interactions with a "x" symbol
#' plot(fit)
#'
#'
#' # Predict using the fitted model
#' ntest=500
#' xtest = matrix(rnorm(ntest*p),ntest,p)
#' ztest =matrix(rnorm(ntest*nz),ntest,nz)
#'
#' pred= predict(fit,xtest,ztest)
#'
#'  #Example where z is not observed in the test set,  but predicted from a  supervised
#'   #  learning method
#'
#'  library(glmnet)
#' n = 20; p = 3 ;nz=3
#' x = matrix(rnorm(n*p), n, p)
#' z =matrix(rnorm(n*nz),n,nz)
#' y = x[,1] +x[,1]*z[,3]+ rnorm(n)
#'fit = pliable(x,z,y)
#' # predict z  from x; here we use glmnet, but any other supervised learning method
#' # could be used
#' zfit=glmnet(x,z,family="mgaussian")
#'
#' # Predict using the fitted model
#' ntest=500
#' xtest = matrix(rnorm(ntest*p),ntest,p)
#' ztest =predict(zfit,xtest,s=zfit$lambda.min)[,,20]
#'
#' pred= predict(fit,xtest,ztest)
#'
#' @export
predict.pliable <- function(object,x,z, type=c("link","response", "coefficients"), lambda=NULL, verbose=FALSE, ...){
       # note- current version uses closest lambda to the values used in the fitted model;
        #  in the future, linear interpolation or exact recomputation should be used
    type = match.arg(type)

    lambda.arg=lambda

    if(is.null(lambda.arg))
      { lambda=object$lambda;isel=1:length(lambda)}

    if(!is.null(lambda.arg)){

          isel=as.numeric(knn1(matrix(object$lambda,ncol=1),matrix(lambda.arg,ncol=1),1:length(object$lambda)))

    }

    if(type=="link" | type=="response"){


if(!is.matrix(z)) z=matrix(z,ncol=1)


    yh=matrix(NA,nrow(x),length(isel))
    ii=0
    for(m in isel){
        ii=ii+1
        if(verbose)  cat(c("m=",m),fill=T)
        if(ncol(object$lp)>0){
            yh[,ii]=predict.pliable1(object,x,z,m,type=type)
            }
            }

return(yh)
}
     if(type=="coefficients"){
        beta=object$beta[,isel]
        theta=object$theta[,,isel]
        betaz=object$betaz[,isel]
        return(list(betaz=betaz,beta=beta,theta=theta))
       }
}
predict.pliable1=function(object,x,z,m,type=c("link","response", "coefficients")){
     type = match.arg(type)
       ni=ncol(x)
       nz=ncol(z)
          n=nrow(x)

   #  cat(colMeans(x),fill=T)


       xc=scale(x,object$xbar,F)
       zc=scale(z,object$zbar,F)



    yhat=object$ybar+object$a0[m]+ zc%*%object$betaz[,m]+xc%*%object$beta[,m]  #note that object$ybar is zero for binomial
  #  yhat=yhat+fit$zbar*fit$betaz[,m]+  sum(fit$xbar*fit$beta[,m])

    ind=which(rowSums(object$theta[,,m,drop=F]!=0)>0)

    for(j in ind){
        xz=matrix(xc[,j],nrow(z),ncol(z))*zc
        yhat=yhat+xz%*%object$theta[j,,m]
       # yhat=yhat+object$xbar[j]*z%*%object$theta[j,,m]
       # yhat=yhat+x[,j]*sum(object$zbar*object$theta[j,,m])
       # yhat=yhat-object$xbar[j]*sum(object$zbar*object$theta[j,,m])
    }
       if(type=="response" & object$family=="binomial") yhat=1/(1+exp(-yhat))
        if(type=="response" & object$family=="cox") yhat=exp(yhat)
    return(yhat)
}

pliableMulti=function(x,z,y,lambda=NULL,nlambda=50,family="binomial",alpha=.5,lambda.min.ratio=NULL,w=rep(1,length(y)),thr=1e-5,tt=NULL, maxit=1000,mxthit=100,mxkbt=100,
                mxth=1000,
                kpmax=100000,kthmax=100000,verbose=TRUE, maxinter=1000, zlinear=TRUE){
    nclass=length(table(y))
    fit=vector("list",nclass)
    for(k in 1:nclass){
        yy=1*(y==k)
        if(k==1) {lambda=lambda; nlambda=nlambda}
        if(k>1) {lambda=fit[[k-1]]$lambda; nlambda=NULL}
        fit[[k]]=pliable(x,z,yy,family=family,lambda=lambda,nlambda=nlambda,alpha=alpha,lambda.min.ratio=lambda.min.ratio,w=w,
                        thr=thr,tt=tt, maxit=maxit,mxthit=mxthit,mxkbt=mxkbt,  mxth=mxth,
                        kpmax=kpmax,kthmax=kthmax,verbose=verbose, maxinter=maxinter, zlinear=zlinear)
    }
    class(fit)="pliableMulti"
    return(fit)
    }



predict.pliableMulti=

    function(fit,x,z, type=c("response", "coefficients"), lambda=NULL, verbose=FALSE){
       # note- current version uses closest lambda to the values used in the fitted model;
        #  in the future, linear interpolation or exact recomputation should be used
    type = match.arg(type)
    nclass=length(fit)
        out=vector("list",nclass)

        for(k in 1:nclass){
            out[[k]]=predict.pliable(fit[[k]],x,z,type=type,lambda=lambda,verbose=verbose)
        }
        if(type=="coefficients") return(out)
        if(type=="response"){
            yhat=array(NA,c(nrow(out[[k]]),ncol(out[[k]]),nclass))
            for(k in 1:nclass){
                yhat[,,k]=out[[k]]
            }
            return(yhat)
            }}







critjerry=function(bb,x,z,y,k,alpha=.5){
    twonorm=function(x){ sqrt(sum(x*x))}
    n=length(y)
    p=ncol(x)
    nz=ncol(z)
    lam=bb$lam[k]
     beta2=bb$beta[,k]
    yhat=z%*%bb$betaz[,k]+x%*%beta2

    theta=bb$theta[,,k]
    if(!is.matrix(theta)) theta=matrix(theta,ncol=1)

    for(j in 1:p){
        xz=matrix(x[,j],nrow(z),ncol(z))*z
        yhat=yhat+xz%*%theta[j,]
    }

    pen=0
    for(j in 1:p){
        pen1=twonorm(theta[j,])+twonorm(c(beta2[j],theta[j,]))
        pen=pen+(1-alpha)*lam*pen1+alpha*lam*sum(abs(theta[j,]))
    }
    out=(1/(2*n))*sum((y-yhat)^2) +pen
return(list(out=out,pen=pen))
}




#' Carries out cross-validation for  a  pliable lasso model over a path of regularization values
#' @param fit  object returned by the pliable function
#' @param x  n by p matrix of predictors
#' @param z n by nz matrix of modifying variables. Note that z may be observed or the predictions from a supervised learning
#'   algorithm that predicts z from validation set features x. IMPORTANT: in the latter case, the z matrix  must be  pre-computed before
#'  calling cv.pliable, using the same CV folds used by cv.pliable. See the example below

#' @param y n-vector of responses. All variables are centered in the function.
#' @param nfolds  number of cross-validation folds
#' @param foldid  vector with values in 1:K, indicating folds for K-fold CV. Default NULL
#' @param keep  Should pre-validated fits be returned?  Default FALSE
#' @param type.measure Error measure for cross-validation: "deviance" (mse) for gaussian family, "deviance" or "class" (misclassification error) for
#'     binomial family
#' @param verbose  Should information  be printed along the way?  Default FALSE
#' @return  predicted values
#' #'
#' @examples \donttest{
#' # Train and cross-validate a pliable lasso model- Gaussian case
#' n = 20 ; p = 3 ;nz=3
#' x = matrix(rnorm(n*p), n, p)
#' mx=colMeans(x)
#' sx=sqrt(apply(x,2,var))
#' x=scale(x,mx,sx)
#' z =matrix(rnorm(n*nz),n,nz)
#' mz=colMeans(z)
#' sz=sqrt(apply(z,2,var))
#' z=scale(z,mz,sz)
#' y =4*x[,1] +5*x[,1]*z[,3]+ 3*rnorm(n)

#fit the model
#' fit = pliable(x,z,y)
#'
#' cvfit=cv.pliable(fit,x,z,y,nfolds=5)
#'   plot(cvfit)
#'
#' # Example of categorical z with 4 levels
#' n = 20; p = 3 ;nz=3
#' x = matrix(rnorm(n*p), n, p)
#'  mx=colMeans(x)
#' sx=sqrt(apply(x,2,var))
#' x=scale(x,mx,sx)
#' z =sample(1:4,size=n,replace=T)
#' zi=model.matrix(~as.factor(z)-1)
#' y = x[,1] +x[,1]*zi[,3]-2*x[,1]*zi[,4]+rnorm(n)
#'fit = pliable(x,zi,y)
#'
#' # Train and cross-validate a pliable lasso model- Binomial case
#' n = 20; p = 3 ;nz=3
#' x = matrix(rnorm(n*p), n, p)
#' mx=colMeans(x)
#' sx=sqrt(apply(x,2,var))
#' x=scale(x,mx,sx)
#' z =matrix(rnorm(n*nz),n,nz)
#' mz=colMeans(z)
#' sz=sqrt(apply(z,2,var))
#' z=scale(z,mz,sz)
#' y =4*x[,1] +5*x[,1]*z[,3]+ 3*rnorm(n)
#' y= 1*(y>0)
#' fit = pliable(x,z,y)
#'
#' cvfit=cv.pliable(fit,x,z,y,type.measure="class", nfolds=5)
#'   plot(cvfit)
#'
#'
#'  # Example where z is not observed in the test set,
#'    #   but predicted from a supervised learning algorithm
#' # NOT RUN
#'
#'#library(glmnet)
#'# n = 20 ; p = 4 ;nz=3
#'# x = matrix(rnorm(n*p), n, p)
#'# mx=colMeans(x)
#'# sx=sqrt(apply(x,2,var))
#'# x=scale(x,mx,sx)
#'# z =matrix(rnorm(n*nz),n,nz)
#'# mz=colMeans(z)
#'# sz=sqrt(apply(z,2,var))
#'# z=scale(z,mz,sz)
#' #y =4*x[,1] +5*x[,1]*z[,3]+ 3*rnorm(n)
#' #fit = pliable(x,z,y)
#'
#'  #nfolds=5
#' #  set.seed(3321)
#' #foldid = sample(rep(seq(nfolds), length = n))
#'
#' # predict z  from x; here we use glmnet, but any other supervised learning procedure
#'   #  could be used
#' #zhat=matrix(NA,n,ncol(z))
#' #for(ii in 1:nfolds){
#' #  zfit=cv.glmnet(x[foldid!=ii,],z[foldid!=ii,],family="mgaussian")
#' #  zhat[foldid==ii,]=predict(zfit,x[foldid==ii,],s=zfit$lambda.min)
#' #}
#'
#' #NOTE that the same foldid vector must be passed to cv.pliable
#' #cvfit=cv.pliable(fit,x,zhat,y,foldid=foldid)
#' #plot(cvfit)
#' #}
#'

#' @export
cv.pliable <-
function(fit,x,z,y,nfolds=10,foldid=NULL,keep=F,type.measure=c("deviance","class"),verbose=TRUE){
    debug=F
    status=NULL # Cox model not yet implemented
    if(debug){

        nfolds=10
        foldid=NULL
        keep=F
        type.measure=c("deviance","class")
        verbose=TRUE
    }
        #put data in time order- need to do this, even though pliable does it too
    if(fit$family=="cox") {
       o=order(y-.0001*status)
        x= x[o,]
        z=z[o,]
        y=y[o]
        status=status[o]

    }

    type.measure=match.arg(type.measure)
   if(fit$family=="gaussian") errfun=errfun.gaussian
    if(fit$family=="binomial" & type.measure=="deviance" )  errfun=errfun.binomial
    if(fit$family=="binomial" & type.measure=="class" )  errfun=function(y,yhat,w) 1*(y!=yhat)*w

  #  if(fit$family=="cox")   errfun=function(yhat,coxinfo) {devc(no,coxinfo$kq, coxinfo$ddq, yhat, coxinfo$iriskq, status)}

    if(fit$family=="cox" & is.null(status)) stop("Must supply status arg with family='cox'")
    BIG=10e9
   ni=ncol(x)
    no=length(y)
    nz=ncol(z)

if(fit$family=="cox") coxinfo.folds=vector("list",nfolds)

    ggg=vector("list",nfolds)

    yhat=array(NA,c(no,length(fit$lambda)))

      if(is.null(foldid))  foldid = sample(rep(seq(nfolds), length = no))
    nfolds=length(table(foldid))

    status.in=NULL
    for(ii in 1:nfolds){
        oo=foldid==ii
     if(fit$family=="cox") status.in=status[!oo]
     if(verbose)  cat(c("\n","FOLD=",ii),fill=T)

      ggg[[ii]]=pliable(x[!oo,,drop=F],z[!oo,,drop=F],y[!oo],family=fit$family,lambda=fit$lambda,alpha=fit$args$alpha ,lambda.min.ratio=fit$args$lambda.min.ratio,w=fit$w[!oo],thr=fit$args$thr,maxit=fit$args$maxit, tt=fit$args$tt, mxthit=fit$args$mxthit ,mxkbt=fit$args$mxkbt, mxth=fit$args$mxth,kpmax=fit$args$kpmax,kthmax=fit$args$kthmax,screen=fit$args$screen)



    if(fit$family=="gaussian" | fit$family=="binomial") yhat[oo,]=predict.pliable(ggg[[ii]],x[oo,,drop=F],z[oo,])
 if(fit$family=="cox") coxinfo.folds[[ii]]=ggg[[ii]]$coxinfo
      }

        nonzero=colSums(fit$beta!=0)

    yhat.preval=NULL

    ym=array(y,dim(yhat))
    if(type.measure=="class") yhat=1*(yhat>mean(y))

    if(fit$family=="gaussian" | fit$family=="binomial") {
         err=errfun(ym,yhat,fit$w)
 cvm=apply(err,2,mean,na.rm=T)
    nn=apply(!is.na(err),2,sum,na.rm=T)
    cvsd=sqrt(apply(err,2,var,na.rm=T)/nn)
    if(keep) yhat.preval=yhat
    }

    if(fit$family=="cox") {
        err=matrix(NA,nfolds,length(fit$lambda))
        for(ii in 1:nfolds){
            oo=foldid==ii
                fit1=predict.pliable(ggg[[ii]],x,z)
                fit2=predict.pliable(ggg[[ii]],x[!oo,,drop=F],z[!oo,])
                for(k in 1:length(fit$lambda)){
                  dev1= devc(no, fit$coxinfo$kq, fit$coxinfo$ddq, fit1[,k], fit$coxinfo$iriskq, status)
                  dev2= devc(sum(!oo), coxinfo.folds[[ii]]$kq, coxinfo.folds[[ii]]$ddq, fit2[,k], coxinfo.folds[[ii]]$iriskq, status[!oo])
                  err[ii,k]=dev1-dev2
                  }

        }
       cvm=apply(err,2,mean,na.rm=T)
          nn=apply(!is.na(err),2,sum,na.rm=T)
    cvsd=sqrt(apply(err,2,var,na.rm=T)/nn)
        }
    cvm.nz=cvm; cvm.nz[nonzero==0]=BIG
    imin=which.min(cvm.nz)


    imin.1se=which(cvm< cvm[imin]+cvsd[imin])[1]

    out=list(lambda=fit$lambda,cvm=cvm,cvsd=cvsd,cvup = cvm +
                                                     cvsd, cvlo = cvm - cvsd, nz=nonzero, df=nonzero,yhat.preval=yhat.preval,lambda.min=fit$lambda[imin],
             lambda.1se=fit$lambda[imin.1se],name="Error")

    class(out)="cv.pliable"
return(out)
}




plot.pliable=
function (x, xvar = c("norm", "lambda", "dev"), label =TRUE,
    ...)
{
    xvar = match.arg(xvar)

    plotCoef(x$beta, theta=x$theta, lambda = x$lambda, df = x$df, dev = x$dev.ratio,
             label = label, xvar = xvar, ...)

}

plot.cv.pliable =
function (x, sign.lambda = 1, ...)
{
    cvobj = x
    xlab = "log(Lambda)"
    if (sign.lambda < 0)
        xlab = paste("-", xlab, sep = "")
    plot.args = list(x = sign.lambda * log(cvobj$lambda), y = cvobj$cvm,
        ylim = range(cvobj$cvup, cvobj$cvlo), xlab = xlab, ylab = cvobj$name,
        type = "n")
    new.args = list(...)
    if (length(new.args))
        plot.args[names(new.args)] = new.args
    do.call("plot", plot.args)
    error.bars(sign.lambda * log(cvobj$lambda), cvobj$cvup, cvobj$cvlo,
               width = 0.01, col = "darkgrey")

    points(sign.lambda * log(cvobj$lambda), cvobj$cvm, pch = 20,
        col = "red")
    axis(side = 3, at = sign.lambda * log(cvobj$lambda), labels = paste(cvobj$nz),
        tick = FALSE, line = 0)
    abline(v = sign.lambda * log(cvobj$lambda.min), lty = 3)
    abline(v = sign.lambda * log(cvobj$lambda.1se), lty = 3)
    invisible()
}
error.bars <-function(x, upper, lower, width = 0.02, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}
plotCoef=
function (beta, theta, norm, lambda, df, dev, label = FALSE, xvar = c("norm","lambda",
"dev"), xlab = iname, ylab = "Coefficients", ...)
{
    which = nonzeroCoef(beta)
    nwhich = length(which)
    switch(nwhich + 1, `0` = {
        warning("No plot produced since all coefficients zero")
        return()
    }, `1` = warning("1 or less nonzero coefficients; glmnet plot is not meaningful"))
    sbeta=beta
    beta = as.matrix(beta[which, , drop = FALSE])
    xvar = match.arg(xvar)
    switch(xvar, norm = {
        index = if (missing(norm)) {apply(abs(beta), 2, sum)+apply(abs(theta),3,sum)} else norm
        iname = "L1 Norm"
        approx.f = 1
    }, lambda = {
        index = log(lambda)
        iname = "Log Lambda"
        approx.f = 0
    }, dev = {
        index = dev
        iname = "Fraction Deviance Explained"
        approx.f = 1
    })
    dotlist = list(...)
    type = dotlist$type
    if (is.null(type))
        matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab,
            type = "l", ...)
    else matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab,
        ...)
    atdf = pretty(index)
    prettydf = approx(x = index, y = df, xout = atdf, rule = 2,
        method = "constant", f = approx.f)$y
    axis(3, at = atdf, labels = prettydf, tcl = NA)
    if (label) {
        nnz = length(which)
        xpos = max(index)
        pos = 4
        if (xvar == "lambda") {
            xpos = min(index)
            pos = 2
        }
        xpos = rep(xpos, nnz)
        ypos = beta[, ncol(beta)]
        text(xpos, ypos, paste(which), cex = 0.5, pos = pos)
    }

     act=which(rowSums(abs(beta))>0)
    ntheta=  apply(theta!=0,c(1,3),sum)

    for(j in act){
        for(i in 1:length(index)){

           if(ntheta[j,i]>0) text(index[i],sbeta[j,i],label="x",cex=.7)
            }}

}

print.pliable=function(x,digits = max(3, getOption("digits") - 3),...){
      cat("\nCall: ", deparse(x$call), "\n\n")

      print(cbind(Lambda=signif(x$lambda,digits), "%Dev"=signif(x$dev.ratio,digits),"Num of main effects"=x$nbeta,"Num with ints"=x$nbeta.with.int,"Tot num of ints"=x$ntheta
                  ))
    }

#' Return estimated coefficients from a fitted pliable model.


#' @param object object returned from a call to pliable
#' @param lambda  values of lambda at whcih predictions are desired. If NULL (default), the path of lambda values from the fitted model.
#'   are used. If lambda is  not NULL,  the predictions are made at the closest values to lambda in the lambda path  from  the fitted model;
#' @param ... Further arguments (not used)
#.
#'
#'  @return  estimated parameters

coef.pliable=
function (object, lambda = NULL,...)
predict(object, lambda=lambda, type = "coefficients")





encode=function(x,vals=NULL){
if(is.null(vals)) vals=names(table(x))
 xx=matrix(0,length(x),length(vals))
 for(i in 1:length(x)){ xx[i,which(x[i]==vals)]=1}
 dimnames(xx)=list(NULL,vals)
 return(xx)
}





nonzeroCoef=
function (beta, bystep = FALSE)
{
    nr = nrow(beta)
    if (nr == 1) {
        if (bystep)
            apply(beta, 2, function(x) if (abs(x) > 0)
                1
            else NULL)
        else {
            if (any(abs(beta) > 0))
                1
            else NULL
        }
    }
    else {
        beta = abs(beta) > 0
        which = seq(nr)
        ones = rep(1, ncol(beta))
        nz = as.vector((beta %*% ones) > 0)
        which = which[nz]
        if (bystep) {
            if (length(which) > 0) {
                beta = as.matrix(beta[which, , drop = FALSE])
                nzel = function(x, which) if (any(x))
                  which[x]
                else NULL
                which = apply(beta, 2, nzel, which)
                if (!is.list(which))
                  which = data.frame(which)
                which
            }
            else {
                dn = dimnames(beta)[[2]]
                which = vector("list", length(dn))
                names(which) = dn
                which
            }
        }
        else which
    }
}

errfun.gaussian=function(y,yhat,w=rep(1,length(y))){  ( w*(y-yhat)^2) }

errfun.binomial=function(y,yhat,w=rep(1,length(y))){
     prob_min = 1e-05
     prob_max = 1 - prob_min
       predmat = pmin(pmax(yhat, prob_min), prob_max)

  -2*w*(y*log(predmat)+(1-y)*log(1-predmat))
}
