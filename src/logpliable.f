c     mortran 2.0     (version of 7/04/75 mod 7/4/87 (ajc))             
      subroutine logpliable(no,ni,nz,x,y,w,alpha,nlam,ulam,kbos,maxinter
     *, linenter,thr,pf,maxit,tbk,mxthit,mxkbt,mxth,kpmx,kthmx,mlam,a0, 
     *lp,istor,fstor,thstor,jerr)
      implicit double precision(a-h,o-z)                                
      double precision y(no),x(no,(ni+nz)),w(no),ulam(nlam),a0(nlam), fs
     *tor(kpmx),thstor(kthmx)
      integer lp(2,nlam),istor(2,kpmx)                                  
      double precision b0(nz+1),fnewb(nz+1),bold(nz+1),gold(nz+1), delta
     *(nz+1),pf(ni)
      double precision, dimension (:,:), allocatable :: theta,z,scrat   
      integer, dimension (:), allocatable :: it                         
      double precision, dimension (:), allocatable :: a,r,rold,rbj,xv   
      double precision, dimension (:), allocatable :: zz,warg,pr        
      double precision, dimension(:), allocatable :: eta,eta0           
      allocate(theta(1:nz,1:mxth),stat=jerr)                            
      if(jerr.ne.0) return                                              
      allocate(z(1:no,1:nz),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      allocate(scrat(1:nz,1:4),stat=jerr)                               
      if(jerr.ne.0) return                                              
      nt=ni+nz                                                          
      allocate(a(1:nt),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(r(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(rold(1:no),stat=jerr)                                    
      if(jerr.ne.0) return                                              
      allocate(rbj(1:no),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(it(1:nt),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xv(1:nt),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(zz(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(warg(1:no),stat=jerr)                                    
      if(jerr.ne.0) return                                              
      allocate(pr(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(eta(1:no),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(eta0(1:no),stat=jerr)                                    
      if(jerr.ne.0) return                                              
      if(kbos .eq. 0)goto 10021                                         
      call dblepr(' number         lambda    var passes  inters', -1,0,0
     *)
      nlp0=0                                                            
10021 continue                                                          
      do 10031 k=1,nz                                                   
      z(:,k)=x(:,ni+k)                                                  
10031 continue                                                          
      continue                                                          
      az=0.0                                                            
      a=0.0                                                             
      it=0                                                              
      nlp=0                                                             
      nth=nlp                                                           
      kp=nth                                                            
      kth=kp                                                            
      r=y                                                               
      sw=sum(w)                                                         
      oma=1.0-alpha                                                     
      do 10041 j=1,nt                                                   
      s=0.0                                                             
      do 10051 i=1,no                                                   
      s=s+w(i)*x(i,j)**2                                                
10051 continue                                                          
      continue                                                          
      xv(j)=s/sw                                                        
10041 continue                                                          
      continue                                                          
      warg=w                                                            
      zz=y                                                              
      eta=0                                                             
      pr=1/(1+exp(-eta))                                                
      w=0.25                                                            
      w=warg*w                                                          
      w=no*w/sum(w)                                                     
      sw=sum(w)                                                         
      a0=sum(w*y)/no                                                    
      mlam=0                                                            
      do 10061 m=1,nlam                                                 
      alm=ulam(m)                                                       
      continue                                                          
10071 continue                                                          
      continue                                                          
10081 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 10091 j=0,nt                                                   
      if(j .ne. 0)goto 10111                                            
      aj=az                                                             
      del=dot_product(w,r)/sw                                           
      az=az+del                                                         
      r=r-del                                                           
      dlx=max(del**2,dlx)                                               
      goto 10091                                                        
10111 continue                                                          
      if(j .le. ni)goto 10131                                           
      if(linenter.eq.0)goto 10091                                       
      aj=a(j)                                                           
      del=dot_product(w*r,x(:,j))/(xv(j)*sw)                            
      a(j)=a(j)+del                                                     
      r=r-del*x(:,j)                                                    
      dlx=max(xv(j)*del**2,dlx)                                         
      goto 10091                                                        
10131 continue                                                          
      almj=alm*pf(j)                                                    
      rbj=r+a(j)*x(:,j)                                                 
      if(it(j) .le. 0)goto 10151                                        
      itj=it(j)                                                         
      do 10161 i=1,no                                                   
      s=0.0                                                             
      do 10171 k=1,nz                                                   
      s=s+theta(k,itj)*x(i,j)*z(i,k)                                    
10171 continue                                                          
      continue                                                          
      rbj(i)=rbj(i)+s                                                   
10161 continue                                                          
      continue                                                          
10151 continue                                                          
      gj=sum(w*rbj*x(:,j))/sw                                           
      s=0.0                                                             
      do 10181 k=1,nz                                                   
      t=0.0                                                             
      do 10191 i=1,no                                                   
      t=t+w(i)*rbj(i)*x(i,j)*z(i,k)                                     
10191 continue                                                          
      continue                                                          
      t=sign(max(0d0,dabs(t/sw)-alpha*almj),t)                           
      s=s+t**2                                                          
10181 continue                                                          
      continue                                                          
      if(sqrt(s) .gt. 2.0*oma*almj*pf(j) .or. abs(gj) .gt. oma*almj)goto
     * 10211
      if(it(j).gt.0) it(j)=-it(j)                                       
      a(j)=0.0                                                          
      goto 10091                                                        
10211 continue                                                          
      aj=0.0                                                            
      do 10221 i=1,no                                                   
      aj=aj+w(i)*x(i,j)*rbj(i)                                          
10221 continue                                                          
      continue                                                          
      aj=aj/sw                                                          
      aj=sign(max(0d0,dabs(aj)-oma*almj),aj)/xv(j)                       
      s=0.0                                                             
      do 10231 k=1,nz                                                   
      t=0.0                                                             
      do 10241 i=1,no                                                   
      t=t+w(i)*(rbj(i)-aj*x(i,j))*x(i,j)*z(i,k)                         
10241 continue                                                          
      continue                                                          
      t=sign(max(0d0,dabs(t/sw)-alpha*almj),t)                           
      s=s+t**2                                                          
10231 continue                                                          
      continue                                                          
      if(sqrt(s) .gt. oma*almj)goto 10261                               
      if(it(j).gt.0) it(j)=-it(j)                                       
      del=aj-a(j)                                                       
      a(j)=aj                                                           
      r=r-del*x(:,j)                                                    
      dlx=max(dlx,xv(j)*del**2)                                         
      goto 10091                                                        
10261 continue                                                          
      if(it(j) .ne. 0)goto 10281                                        
      nth=nth+1                                                         
      if(nth .le. mxth)goto 10301                                       
      jerr=100                                                          
      return                                                            
10301 continue                                                          
      it(j)=nth                                                         
      b0(2:(nz+1))=0.0                                                  
      goto 10271                                                        
10281 if(it(j) .ge. 0)goto 10311                                        
      it(j)=-it(j)                                                      
      b0(2:(nz+1))=0.0                                                  
      goto 10321                                                        
10311 continue                                                          
      b0(2:(nz+1))=theta(:,it(j))                                       
10321 continue                                                          
10271 continue                                                          
      b0(1)=a(j)                                                        
      bold=b0                                                           
      rssold=sum(w*r**2)/sw                                             
      rold=r                                                            
      kerr=0                                                            
      do 10331 kthit=1,mxthit                                           
      grbeta=0.0                                                        
      do 10341 i=1,no                                                   
      grbeta=grbeta-w(i)*r(i)*x(i,j)                                    
10341 continue                                                          
      continue                                                          
      grbeta=grbeta/sw                                                  
      do 10351 k=1,nz                                                   
      s=0.0                                                             
      do 10361 i=1,no                                                   
      s=s-w(i)*r(i)*x(i,j)*z(i,k)                                       
10361 continue                                                          
      continue                                                          
      scrat(k,4)=s/sw                                                   
10351 continue                                                          
      continue                                                          
      if(kthit .ne. 1)goto 10381                                        
      gold(1)=grbeta                                                    
      do 10391 k=1,nz                                                   
      gold(k+1)=scrat(k,4)                                              
10391 continue                                                          
      continue                                                          
10381 continue                                                          
      tt=tbk                                                            
      do 10401 kbt=1,mxkbt                                              
      call solveab(nz,b0,grbeta,scrat(:,4),alpha,almj,tt,  scrat(:,1),sc
     *rat(:,2),fnewb,irrflag)
      if(irrflag .eq. 0)goto 10421                                      
      kerr=1                                                            
      goto 10402                                                        
10421 continue                                                          
      r=r-(fnewb(1)-b0(1))*x(:,j)                                       
      do 10431 k=1,nz                                                   
      r=r-(fnewb(k+1)-b0(k+1))*x(:,j)*z(:,k)                            
10431 continue                                                          
      continue                                                          
      rssnew=sum(w*r**2)/sw                                             
      delta=fnewb-bold                                                  
      cri=rssnew-rssold  -2.0*sum(delta*gold)-sum(delta**2)/tt          
      if(cri.le.0.0)goto 10402                                          
      tt=0.9*tt                                                         
10401 continue                                                          
10402 continue                                                          
      if(kerr.eq.1)goto 10332                                           
      if(konv(nz+1,b0,fnewb).ne.0)goto 10332                            
      b0=fnewb                                                          
10331 continue                                                          
10332 continue                                                          
      if(kerr .ne. 1)goto 10451                                         
      if(it(j).gt.0) it(j)=-it(j)                                       
      r=rold                                                            
      goto 10091                                                        
10451 continue                                                          
      aj=a(j)                                                           
      a(j)=fnewb(1)                                                     
      theta(:,it(j))=fnewb(2:(nz+1))                                    
      dlx=max(dlx,xv(j)*(a(j)-aj)**2)                                   
10091 continue                                                          
      continue                                                          
      if(dlx.lt.thr)goto 10082                                          
      if(nlp .le. maxit)goto 10471                                      
      jerr=-m                                                           
      return                                                            
10471 continue                                                          
      goto 10081                                                        
10082 continue                                                          
      eta0=eta                                                          
      eta=zz-r                                                          
      del2=sum(abs(eta-eta0))/no                                        
      pr=1/(1+exp(-eta))                                                
      zz=eta+(y-pr)/w                                                   
      w=warg*w                                                          
      w=no*w/sum(w)                                                     
      sw=sum(w)                                                         
      a0(m)=az                                                          
      r=zz-eta                                                          
      if(del2.lt.0.01)goto 10072                                        
      goto 10071                                                        
10072 continue                                                          
      a0(m)=az                                                          
      lp(1,m)=kp+1                                                      
      do 10481 j=1,nt                                                   
      if(a(j).eq.0.0)goto 10481                                         
      call modstor(j,a(j),it(j),nz,theta(:,it(j)),kpmx,kthmx,  istor,fst
     *or,thstor,kp,kth,jerr)
      if(jerr .eq. 0)goto 10501                                         
      jerr=200                                                          
      return                                                            
10501 continue                                                          
10481 continue                                                          
      continue                                                          
      lp(2,m)=kp                                                        
      jiter=0                                                           
      do 10511 j=1,ni                                                   
      if(it(j).le.0)goto 10511                                          
      do 10521 k=1,nz                                                   
      if(theta(k,it(j)).ne.0.0) jiter=jiter+1                           
10521 continue                                                          
      continue                                                          
10511 continue                                                          
      continue                                                          
      if(kbos .eq. 0)goto 10541                                         
      call dblepr(' ',-1,0,0)                                           
      call intpr('Step=',-1,m,1)                                        
      call dblepr('Lambda=',-1,ulam(m),1)                               
      call intpr('Number of intns=',-1,jiter,1)                         
      call dblepr(' ',-1,0,0)                                           
      nlp0=nlp                                                          
10541 continue                                                          
      if(jiter .lt. maxinter)goto 10561                                 
      if(kbos .eq. 0)goto 10581                                         
      call dblepr('Interaction count exceeded',-1,0,0)                  
10581 continue                                                          
      mlam=mlam+1                                                       
      goto 10062                                                        
10561 continue                                                          
      mlam=mlam+1                                                       
10061 continue                                                          
10062 continue                                                          
      deallocate(theta,z,it,a,r,rold,rbj,scrat,xv,zz,warg,pr,eta,eta0)  
      return                                                            
      end                                                               
