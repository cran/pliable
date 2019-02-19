"
c
c                          pliable (2/27/18)
c                            
c                        
c                 Pliable lasso squared-error loss
c
c
c call pliable(no,ni,nz,x,y,w,alpha,nlam,ulam,kbos,maxinter,linenter,
c   thr,pf,maxit,tbk,mxthit,mxkbt,mxth,kpmx,kthmx,mlam,a0,lp,istor,fstor,
c   thstor,jerr)
c
c input:
c
c   no = number of observations
c   ni = number of predictor variables
c   nz = number of interacting variables
c   x(no,ni+nz) = predictor and interacting variables data matrix
c      x(no,1:ni) = predictor variables
c      x(no,(ni+1):(ni+nz)) = interacting variables
c     (all columns centered)
c   y(no) = response vector (centered)
c   w(no)= observation weights
c   alpha = interaction strength parameter (in [0,1])
c   nlam = number of lamda values
c   ulam(nlam) = user supplied lamda values in descending order
c   kbos = 0/1 => don't/do print output while executing
c   maxinter = maximum number of interacting variables allowed
c   linenter = 0/1 => don't/do incorporate linear terms for interacting variables
c   thr = convergence threshold for each lamda solution.
c      (suggested value, thr=1.0e-5)
C   pf: penalty factors for each X variable. Default: vector of 1s
c   maxit = maximum allowed number of passes over the data for all lambda
c      values (suggested value, maxit = 100000)
c   tbk = back tracking parameter for theta solution (suggested value, tbk=1.0)
c   mxthit = maximum number of theta solution iterations (suggested value, mxthit=100)
c   mxkbt = maximum number of backtracking iterations (suggested value, mxthit=100)
c   mxth = maximum internal theta storage
c   kpmax = maximum dimension of istor,rstor (see below)
c   kthmax = maximum dimension of thstor (see below)
c
c output:
c
c mlam = number of lambda solutions returned
c a0(nlam) = intercept for each lambda solution
c lp(2,nlam) = pointers to model for each lambda value
c istor(2,kpmx),fstor(kpmx),thstor(kthmax) = model storage for all solutions
c jerr = 0/1 => no/yes storage error
c 
c
c                            Utility Routines
c
c                  Model predictions for kth lambda value solution
c
c call modpred(no,ni,nz,x,a0(k),lp(:,k),istor,fstor,thstor,fh)

c input:
c
c   no = number of observations to be predicted
c   ni = number of x-variables
c   nz = number of z-variables
c   x(no,ni+nz) = x & z input data matrix
c   a0(k) = intercept
c   lp,istor,fstor,thstor = output from pliable
c   
c output
c
c   fh(no) = observation model predictions
c
c 
c                   Coefficients for kth lambda solution
c
c call modsoln(nz,lp(:,k),istor,fstor,thstor,kv,iv,av,it,kz,th)
c
c input:
c
c   nz = number of z-variables
c   lp,istor,fstor,thstor, = output from pliable
c
c output:
c
c   kv = number of non-zero coefficients
c   iv(kv) = identites of non-zero coefficients
c   av(kv) = corresponding non-zero values
c   it(kv) = interaction pointer for each non-zero coefficient:
c      it(j)=0 => no interaction for iv(j) coefficient
c      it(j)>0 => interaction coefficients (thetas) in th(1:nz,it(j))
c   kz = number of interacting coefficient sets for this solution
c   th(nz,kz) = all interaction coefficient sets for this solution
c
c
c  
"
subroutine pliable(no,ni,nz,x,y,w,alpha,nlam,ulam,kbos,maxinter,
linenter,thr,pf,maxit,tbk,mxthit,mxkbt,mxth,kpmx,kthmx,mlam,a0,
lp,istor,fstor,thstor,jerr);
real y(no),x(no,(ni+nz)),w(no),ulam(nlam),a0(nlam),
fstor(kpmx),thstor(kthmx);
integer lp(2,nlam),istor(2,kpmx);
real b0(nz+1),fnewb(nz+1),bold(nz+1),gold(nz+1),
  delta(nz+1),pf(ni);
%fortran
      real, dimension (:,:), allocatable :: theta,z,scrat
      integer, dimension (:), allocatable :: it
      real, dimension (:), allocatable :: a,r,rold,rbj,xv
      allocate(theta(1:nz,1:mxth),stat=jerr)
      if(jerr.ne.0) return
      allocate(z(1:no,1:nz),stat=jerr)
      if(jerr.ne.0) return
      allocate(scrat(1:nz,1:4),stat=jerr)
      if(jerr.ne.0) return
%mortran
nt=ni+nz;
allocate(a(1:nt),stat=ierr); jerr=jerr+ierr;
allocate(r(1:no),stat=ierr); jerr=jerr+ierr;
allocate(rold(1:no),stat=ierr); jerr=jerr+ierr;
allocate(rbj(1:no),stat=ierr); jerr=jerr+ierr;
allocate(it(1:nt),stat=ierr); jerr=jerr+ierr;
allocate(xv(1:nt),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
if kbos.ne.0 < 
"  <w>; (' number         lambda    var passes  inters');
"
    nlp0=0;>
<k=1,nz; z(:,k)=x(:,ni+k);>
az=0.0; a=0.0; it=0; /nlp,nth,kp,kth/=0; r=y; sw=sum(w); oma=1.0-alpha;
<j=1,nt; s=0.0; <i=1,no; s=s+w(i)*x(i,j)**2;> xv(j)=s/sw;>; mlam=0;
<m=1,nlam;  alm=ulam(m);
   loop < nlp=nlp+1; dlx=0.0;
      <j=0,nt;
         if j.eq.0 < aj=az; del=dot_product(w,r)/sw;
            az=az+del; r=r-del; dlx=max(del**2,dlx);
            next;
         >
         if j.gt.ni < if(linenter.eq.0) next;
            aj=a(j);
            del=dot_product(w*r,x(:,j))/(xv(j)*sw); a(j)=a(j)+del;
            r=r-del*x(:,j); dlx=max(xv(j)*del**2,dlx);
            next;
         >
         almj=alm*pf(j);
         rbj=r+a(j)*x(:,j);
         if it(j).gt.0 <itj=it(j);
            <i=1,no; s=0.0; 
               <k=1,nz; s=s+theta(k,itj)*x(i,j)*z(i,k);>
                  rbj(i)=rbj(i)+s;
            > 
         >  
         gj=sum(w*rbj*x(:,j))/sw;
         s=0.0;
         <k=1,nz; t=0.0; 
            <i=1,no; t=t+w(i)*rbj(i)*x(i,j)*z(i,k);>
            t=sign(max(0d0,abs(t/sw)-alpha*almj),t);
            s=s+t**2;
         >
        "check if both are zero"
         if sqrt(s).le.2.0*oma*almj.and.abs(gj).le.oma*almj <
            if(it(j).gt.0) it(j)=-it(j); a(j)=0.0; next;
         >
         aj=0.0;
         <i=1,no; aj=aj+w(i)*x(i,j)*rbj(i);> aj=aj/sw;
         aj=sign(max(0d0,abs(aj)-oma*almj),aj)/xv(j);
         s=0.0;
         <k=1,nz; t=0.0; 
            <i=1,no; t=t+w(i)*(rbj(i)-aj*x(i,j))*x(i,j)*z(i,k);>
            t=sign(max(0d0,abs(t/sw)-alpha*almj),t);
            s=s+t**2;
         >
         "check if theta =0, with beta ne 0"      
         if sqrt(s).le.oma*almj <
            if(it(j).gt.0) it(j)=-it(j); del=aj-a(j); a(j)=aj;
            r=r-del*x(:,j); dlx=max(dlx,xv(j)*del**2); next;
         >
         "both a(j) and theta(:,j) not zero"
         if it(j).eq.0 < nth=nth+1; if nth.gt.mxth < jerr=100; return;>
            it(j)=nth; b0(2:(nz+1))=0.0;
         >
         elseif it(j).lt.0 < it(j)=-it(j); b0(2:(nz+1))=0.0;>
         else < b0(2:(nz+1))=theta(:,it(j));>
         b0(1)=a(j); bold=b0; rssold=sum(w*r**2)/sw; rold=r; kerr=0;
         <kthit=1,mxthit;
            grbeta=0.0; <i=1,no; grbeta=grbeta-w(i)*r(i)*x(i,j);>
            grbeta=grbeta/sw;
            <k=1,nz; s=0.0;
               <i=1,no; s=s-w(i)*r(i)*x(i,j)*z(i,k);> scrat(k,4)=s/sw;
            >
            if kthit.eq.1 < gold(1)=grbeta; <k=1,nz; gold(k+1)=scrat(k,4);>>
            tt=tbk;
            <kbt=1,mxkbt;
               call solveab(nz,b0,grbeta,scrat(:,4),alpha,almj,tt,
                  scrat(:,1),scrat(:,2),fnewb,irrflag);
               if irrflag.ne.0 < kerr=1; exit;>
               r=r-(fnewb(1)-b0(1))*x(:,j);
               <k=1,nz; r=r-(fnewb(k+1)-b0(k+1))*x(:,j)*z(:,k);>
               rssnew=sum(w*r**2)/sw; delta=fnewb-bold;
               cri=rssnew-rssold
                  -2.0*sum(delta*gold)-sum(delta**2)/tt;
               if(cri.le.0.0) exit; tt=0.9*tt;
            >
            if(kerr.eq.1) exit;
            if(konv(nz+1,b0,fnewb).ne.0) exit;
            b0=fnewb;
         >
         if kerr.eq.1 < if(it(j).gt.0) it(j)=-it(j); r=rold; next;>
         aj=a(j); a(j)=fnewb(1); theta(:,it(j))=fnewb(2:(nz+1));
         dlx=max(dlx,xv(j)*(a(j)-aj)**2);
      > " end of  j=1,nt loop"
      if(dlx.lt.thr) exit;
      if nlp.gt.maxit < jerr=-m; return;>
   >"end of  loop statement"   
   a0(m)=az; lp(1,m)=kp+1;
   <j=1,nt; if(a(j).eq.0.0) next;
      call modstor(j,a(j),it(j),nz,theta(:,it(j)),kpmx,kthmx,
         istor,fstor,thstor,kp,kth,jerr);
      if jerr.ne.0 < jerr=200; return;>
   >
   lp(2,m)=kp;
   jiter=0;
   <j=1,ni; if(it(j).le.0) next;
      <k=1,nz; if(theta(k,it(j)).ne.0.0) jiter=jiter+1;>
   >
   if kbos.ne.0 <
      
      
"      <w> m,ulam(m),nlp-nlp0,jiter; (i4,'        ',g12.4,'  ',i6,'   ',i6);
"
       call dblepr(' ',-1,0,0);
call intpr('Step=',-1,m,1);
      call dblepr('Lambda=',-1,ulam(m),1);
      call intpr('Number of intns=',-1,jiter,1);
        call dblepr(' ',-1,0,0);
nlp0=nlp;
   >  
   if jiter.ge.maxinter <
      if kbos.ne.0 <
"          <w>; ('interaction count exceeded');
"  
            call dblepr('Interaction count exceeded',-1,0,0);
>
      mlam=mlam+1; exit;
   >
   mlam=mlam+1;
>
deallocate(theta,z,it,a,r,rold,rbj,scrat,xv);
return;
end;
function konv(n,u,v);
parameter(eps=1.0e-7);
real u(n),v(n);
scl=0.5*eps*sum(abs(u)+abs(v))/n;
<k=1,n; if abs(u(k)-v(k)).gt.eps < konv=0; return;>>
konv=1;
return;
end;
subroutine modstor(iv,av,noth,nz,thv,kpmx,kthmx,istor,fstor,thstor,kp,kth,
jerr);
real thv(nz),fstor(kpmx),thstor(nz,kthmx); integer istor(2,kpmx);
kp=kp+1; jerr=0; if kp.gt.kpmx < jerr=1; return;>
istor(1,kp)=iv; fstor(kp)=av;
if noth.le.0 < istor(2,kp)=0; return;>
kth=kth+1; if kth.gt.kthmx < jerr=2; return;>
istor(2,kp)=kth; thstor(:,kth)=thv;
return;
end;
subroutine modpred(no,ni,nz,x,a0,kpl,istor,fstor,thstor,fh);
real x(no,ni+nz),fstor(*),thstor(nz,*),fh(no);
integer istor(2,*),kpl(2);
<i=1,no; f=a0;
   <k=kpl(1),kpl(2); is1=istor(1,k); is2=istor(2,k);
      f=f+fstor(k)*x(i,is1);
      if(is2.eq.0) next; xs=x(i,is1);
      <l=1,nz; f=f+thstor(l,is2)*xs*x(i,l+nz);>
   >
   fh(i)=f;
>
return;
end;
subroutine modsoln(nz,kpl,istor,fstor,thstor,kv,iv,av,it,kz,tv);
real fstor(*), thstor(nz,*),av(*),tv(nz,*);
integer kpl(2),istor(2,*),iv(*),it(*);
kv=0; kz=0;
<k=kpl(1),kpl(2); kv=kv+1; iv(kv)=istor(1,k); av(kv)=fstor(k);
   if istor(2,k).eq.0 < it(kv)=0; next;>
   kz=kz+1; it(kv)=kz; is2=istor(2,k);
   <l=1,nz; tv(l,kz)=thstor(l,is2);>
>
return;
end;
%%
