 #helper functions. No need to export

calcz=function(x,kq,irindq,icensq,iriskq,ddq,beta){
n<-nrow(x)
etahat<-x%*%beta

#load.if.needed("/u/tibs/PAPERS/lasso/cox/cox.o")
iriskq<-c(iriskq,rep(0,n-length(iriskq)))
ddq<-c(ddq,rep(0,n-length(ddq)))

mode(etahat)<-"double"
mode(iriskq)<-"integer"
mode(icensq)<-"integer"
mode(ddq)<-"integer"

junk<-.Fortran("calcz",
z=double(n),
wz=double(n),
as.integer(n),
as.integer(kq),
ddq,
etahat,
irindq=integer(n),
iriskq,
icensq,
scrt1=double(n),
scrt2=double(n),
scrt3=double(n),
scrt4=double(n),
scrt5=double(n))

return(list(z=junk$z,wz=junk$wz,grad=junk$scrt4))
}
devc=
function (n, kq, ddq, fits, iriskq, icensq) 
{
    mode(fits) <- "double"
    mode(ddq) <- "integer"
    mode(iriskq) <- "integer"
    mode(icensq) <- "integer"
    junk <- .Fortran("devc", as.integer(n), as.integer(kq), ddq, 
        fits, iriskq, icensq, sfitsq = double(n), scrt1 = double(n), 
        value = double(1))
    return(junk$value)
}
quantitative.func <- function(x, y, s0 = 0) {
        # regression of x on y
        my = mean(y)
        yy <- y - my
        temp <- x %*% yy
        mx = rowMeans(x)
        syy = sum(yy^2)
        scor <- temp/syy
        b0hat <- mx - scor * my
        ym = matrix(y, nrow = nrow(x), ncol = ncol(x), byrow = T)
        xhat <- matrix(b0hat, nrow = nrow(x), ncol = ncol(x)) + ym *
                matrix(scor, nrow = nrow(x), ncol = ncol(x))
        sigma <- sqrt(rowSums((x - xhat)^2)/(ncol(xhat) - 2))
        sd <- sigma/sqrt(syy)
        tt <- scor/(sd + s0)
        return(list(tt = tt, numer = scor, sd = sd))
}


# these are double prec versions of R cox functions, for use in pliable


initcx=function(x,y,ic){
n<-length(y)
kq<-0
px<-ncol(x)


etahat<-rep(0,n)
muhat<-etahat
iriskq<-rep(0,n)
ddq<-rep(0,n)
tq<-rep(0,n)


#dyn.load("/Users/tibs/dropbox/PAPERS/lasso/cox/S.so")
#if(!is.loaded("initcx"))  dyn.load("/Users/tibs/dropbox/PAPERS/weightree2/jerry/pliable2.so")
if(!is.loaded("initcx"))  dyn.load("/Users/tibs/dropbox/PAPERS/weightree2/jerry/cox/pliable2.so")
junk<-.Fortran("initcx",
as.integer(n),
as.integer(px),
as.double(y),
kq=as.integer(kq),
as.integer(ic),
iriskq=as.integer(iriskq),
ddq=as.integer(ddq),
tq=as.double(tq))
kq<-junk$kq;
tq<-junk$tq[1:kq];ddq<-junk$ddq[1:kq];iriskq<-junk$iriskq[1:kq]


return(list(kq=kq,tq=tq,ddq=ddq,iriskq=iriskq))
    }
    quantitative.func <- function(x, y, s0 = 0) {
        # regression of x on y
        my = mean(y)
        yy <- y - my
        temp <- x %*% yy
        mx = rowMeans(x)
        syy = sum(yy^2)
        scor <- temp/syy
        b0hat <- mx - scor * my
        ym = matrix(y, nrow = nrow(x), ncol = ncol(x), byrow = T)
        xhat <- matrix(b0hat, nrow = nrow(x), ncol = ncol(x)) + ym *
                matrix(scor, nrow = nrow(x), ncol = ncol(x))
        sigma <- sqrt(rowSums((x - xhat)^2)/(ncol(xhat) - 2))
        sd <- sigma/sqrt(syy)
        tt <- scor/(sd + s0)
        return(list(tt = tt, numer = scor, sd = sd))
}


# these are double prec versions of R cox functions, for use in pliable


initcx=function(x,y,ic){
n<-length(y)
kq<-0
px<-ncol(x)


etahat<-rep(0,n)
muhat<-etahat
iriskq<-rep(0,n)
ddq<-rep(0,n)
tq<-rep(0,n)


#dyn.load("/Users/tibs/dropbox/PAPERS/lasso/cox/S.so")
#if(!is.loaded("initcx"))  dyn.load("/Users/tibs/dropbox/PAPERS/weightree2/jerry/pliable2.so")
if(!is.loaded("initcx"))  dyn.load("/Users/tibs/dropbox/PAPERS/weightree2/jerry/cox/pliable2.so")
junk<-.Fortran("initcx",
as.integer(n),
as.integer(px),
as.double(y),
kq=as.integer(kq),
as.integer(ic),
iriskq=as.integer(iriskq),
ddq=as.integer(ddq),
tq=as.double(tq))
kq<-junk$kq;
tq<-junk$tq[1:kq];ddq<-junk$ddq[1:kq];iriskq<-junk$iriskq[1:kq]


return(list(kq=kq,tq=tq,ddq=ddq,iriskq=iriskq))
    }
    



    varr <- function(x, meanx = NULL) {
        n <- ncol(x)
        p <- nrow(x)
        Y <- matrix(1, nrow = n, ncol = 1)
        if (is.null(meanx)) {
                meanx <- rowMeans(x)
        }
        ans <- rep(1, p)
        xdif <- x - meanx %*% t(Y)
        ans <- (xdif^2) %*% rep(1/(n - 1), n)
        ans <- drop(ans)
        return(ans)
}
ttest.func <- function(x, y, s0 = 0, sd = NULL) {
        n1 <- sum(y == 1)
        n2 <- sum(y == 2)
        p <- nrow(x)
        m1 <- rowMeans(x[, y == 1, drop = F])
        m2 <- rowMeans(x[, y == 2, drop = F])
        if (is.null(sd)) {
                sd <- sqrt(((n2 - 1) * varr(x[, y == 2], meanx = m2) +
                        (n1 - 1) * varr(x[, y == 1], meanx = m1)) * (1/n1 +
                        1/n2)/(n1 + n2 - 2))
        }
        numer <- m2 - m1
        dif.obs <- (numer)/(sd + s0)
        return(list(tt = dif.obs, numer = numer, sd = sd))
}

