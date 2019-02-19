c     mortran 2.0     (version of 7/04/75 mod 7/4/87 (ajc))             
      subroutine pliable(no,ni,nz,x,y,w,alpha,nlam,ulam,kbos,maxinter, l
     *inenter,thr,pf,maxit,tbk,mxthit,mxkbt,mxth,kpmx,kthmx,mlam,a0, lp,
     *istor,fstor,thstor,jerr)
      implicit double precision(a-h,o-z)                                
      double precision y(no),x(no,(ni+nz)),w(no),ulam(nlam),a0(nlam), fs
     *tor(kpmx),thstor(kthmx)
      integer lp(2,nlam),istor(2,kpmx)                                  
      double precision b0(nz+1),fnewb(nz+1),bold(nz+1),gold(nz+1),  delt
     *a(nz+1),pf(ni)
      double precision, dimension (:,:), allocatable :: theta,z,scrat   
      integer, dimension (:), allocatable :: it                         
      double precision, dimension (:), allocatable :: a,r,rold,rbj,xv   
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
      if(kbos .eq. 0)goto 10021                                         
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
      mlam=0                                                            
      do 10061 m=1,nlam                                                 
      alm=ulam(m)                                                       
      continue                                                          
10071 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 10081 j=0,nt                                                   
      if(j .ne. 0)goto 10101                                            
      aj=az                                                             
      del=dot_product(w,r)/sw                                           
      az=az+del                                                         
      r=r-del                                                           
      dlx=max(del**2,dlx)                                               
      goto 10081                                                        
10101 continue                                                          
      if(j .le. ni)goto 10121                                           
      if(linenter.eq.0)goto 10081                                       
      aj=a(j)                                                           
      del=dot_product(w*r,x(:,j))/(xv(j)*sw)                            
      a(j)=a(j)+del                                                     
      r=r-del*x(:,j)                                                    
      dlx=max(xv(j)*del**2,dlx)                                         
      goto 10081                                                        
10121 continue                                                          
      almj=alm*pf(j)                                                    
      rbj=r+a(j)*x(:,j)                                                 
      if(it(j) .le. 0)goto 10141                                        
      itj=it(j)                                                         
      do 10151 i=1,no                                                   
      s=0.0                                                             
      do 10161 k=1,nz                                                   
      s=s+theta(k,itj)*x(i,j)*z(i,k)                                    
10161 continue                                                          
      continue                                                          
      rbj(i)=rbj(i)+s                                                   
10151 continue                                                          
      continue                                                          
10141 continue                                                          
      gj=sum(w*rbj*x(:,j))/sw                                           
      s=0.0                                                             
      do 10171 k=1,nz                                                   
      t=0.0                                                             
      do 10181 i=1,no                                                   
      t=t+w(i)*rbj(i)*x(i,j)*z(i,k)                                     
10181 continue                                                          
      continue                                                          
      t=sign(max(0d0,abs(t/sw)-alpha*almj),t)                           
      s=s+t**2                                                          
10171 continue                                                          
      continue                                                          
      if(sqrt(s) .gt. 2.0*oma*almj .or. abs(gj) .gt. oma*almj)goto 10201
      if(it(j).gt.0) it(j)=-it(j)                                       
      a(j)=0.0                                                          
      goto 10081                                                        
10201 continue                                                          
      aj=0.0                                                            
      do 10211 i=1,no                                                   
      aj=aj+w(i)*x(i,j)*rbj(i)                                          
10211 continue                                                          
      continue                                                          
      aj=aj/sw                                                          
      aj=sign(max(0d0,abs(aj)-oma*almj),aj)/xv(j)                       
      s=0.0                                                             
      do 10221 k=1,nz                                                   
      t=0.0                                                             
      do 10231 i=1,no                                                   
      t=t+w(i)*(rbj(i)-aj*x(i,j))*x(i,j)*z(i,k)                         
10231 continue                                                          
      continue                                                          
      t=sign(max(0d0,abs(t/sw)-alpha*almj),t)                           
      s=s+t**2                                                          
10221 continue                                                          
      continue                                                          
      if(sqrt(s) .gt. oma*almj)goto 10251                               
      if(it(j).gt.0) it(j)=-it(j)                                       
      del=aj-a(j)                                                       
      a(j)=aj                                                           
      r=r-del*x(:,j)                                                    
      dlx=max(dlx,xv(j)*del**2)                                         
      goto 10081                                                        
10251 continue                                                          
      if(it(j) .ne. 0)goto 10271                                        
      nth=nth+1                                                         
      if(nth .le. mxth)goto 10291                                       
      jerr=100                                                          
      return                                                            
10291 continue                                                          
      it(j)=nth                                                         
      b0(2:(nz+1))=0.0                                                  
      goto 10261                                                        
10271 if(it(j) .ge. 0)goto 10301                                        
      it(j)=-it(j)                                                      
      b0(2:(nz+1))=0.0                                                  
      goto 10311                                                        
10301 continue                                                          
      b0(2:(nz+1))=theta(:,it(j))                                       
10311 continue                                                          
10261 continue                                                          
      b0(1)=a(j)                                                        
      bold=b0                                                           
      rssold=sum(w*r**2)/sw                                             
      rold=r                                                            
      kerr=0                                                            
      do 10321 kthit=1,mxthit                                           
      grbeta=0.0                                                        
      do 10331 i=1,no                                                   
      grbeta=grbeta-w(i)*r(i)*x(i,j)                                    
10331 continue                                                          
      continue                                                          
      grbeta=grbeta/sw                                                  
      do 10341 k=1,nz                                                   
      s=0.0                                                             
      do 10351 i=1,no                                                   
      s=s-w(i)*r(i)*x(i,j)*z(i,k)                                       
10351 continue                                                          
      continue                                                          
      scrat(k,4)=s/sw                                                   
10341 continue                                                          
      continue                                                          
      if(kthit .ne. 1)goto 10371                                        
      gold(1)=grbeta                                                    
      do 10381 k=1,nz                                                   
      gold(k+1)=scrat(k,4)                                              
10381 continue                                                          
      continue                                                          
10371 continue                                                          
      tt=tbk                                                            
      do 10391 kbt=1,mxkbt                                              
      call solveab(nz,b0,grbeta,scrat(:,4),alpha,almj,tt,  scrat(:,1),sc
     *rat(:,2),fnewb,irrflag)
      if(irrflag .eq. 0)goto 10411                                      
      kerr=1                                                            
      goto 10392                                                        
10411 continue                                                          
      r=r-(fnewb(1)-b0(1))*x(:,j)                                       
      do 10421 k=1,nz                                                   
      r=r-(fnewb(k+1)-b0(k+1))*x(:,j)*z(:,k)                            
10421 continue                                                          
      continue                                                          
      rssnew=sum(w*r**2)/sw                                             
      delta=fnewb-bold                                                  
      cri=rssnew-rssold  -2.0*sum(delta*gold)-sum(delta**2)/tt          
      if(cri.le.0.0)goto 10392                                          
      tt=0.9*tt                                                         
10391 continue                                                          
10392 continue                                                          
      if(kerr.eq.1)goto 10322                                           
      if(konv(nz+1,b0,fnewb).ne.0)goto 10322                            
      b0=fnewb                                                          
10321 continue                                                          
10322 continue                                                          
      if(kerr .ne. 1)goto 10441                                         
      if(it(j).gt.0) it(j)=-it(j)                                       
      r=rold                                                            
      goto 10081                                                        
10441 continue                                                          
      aj=a(j)                                                           
      a(j)=fnewb(1)                                                     
      theta(:,it(j))=fnewb(2:(nz+1))                                    
      dlx=max(dlx,xv(j)*(a(j)-aj)**2)                                   
10081 continue                                                          
      continue                                                          
      if(dlx.lt.thr)goto 10072                                          
      if(nlp .le. maxit)goto 10461                                      
      jerr=-m                                                           
      return                                                            
10461 continue                                                          
      goto 10071                                                        
10072 continue                                                          
      a0(m)=az                                                          
      lp(1,m)=kp+1                                                      
      do 10471 j=1,nt                                                   
      if(a(j).eq.0.0)goto 10471                                         
      call modstor(j,a(j),it(j),nz,theta(:,it(j)),kpmx,kthmx,  istor,fst
     *or,thstor,kp,kth,jerr)
      if(jerr .eq. 0)goto 10491                                         
      jerr=200                                                          
      return                                                            
10491 continue                                                          
10471 continue                                                          
      continue                                                          
      lp(2,m)=kp                                                        
      jiter=0                                                           
      do 10501 j=1,ni                                                   
      if(it(j).le.0)goto 10501                                          
      do 10511 k=1,nz                                                   
      if(theta(k,it(j)).ne.0.0) jiter=jiter+1                           
10511 continue                                                          
      continue                                                          
10501 continue                                                          
      continue                                                          
      if(kbos .eq. 0)goto 10531                                         
      call dblepr(' ',-1,0,0)                                           
      call intpr('Step=',-1,m,1)                                        
      call dblepr('Lambda=',-1,ulam(m),1)                               
      call intpr('Number of intns=',-1,jiter,1)                         
      call dblepr(' ',-1,0,0)                                           
      nlp0=nlp                                                          
10531 continue                                                          
      if(jiter .lt. maxinter)goto 10551                                 
      if(kbos .eq. 0)goto 10571                                         
      call dblepr('Interaction count exceeded',-1,0,0)                  
10571 continue                                                          
      mlam=mlam+1                                                       
      goto 10062                                                        
10551 continue                                                          
      mlam=mlam+1                                                       
10061 continue                                                          
10062 continue                                                          
      deallocate(theta,z,it,a,r,rold,rbj,scrat,xv)                      
      return                                                            
      end                                                               
      function konv(n,u,v)                                              
      implicit double precision(a-h,o-z)                                
      parameter(eps=1.0e-7)                                             
      double precision u(n),v(n)                                        
      scl=0.5*eps*sum(abs(u)+abs(v))/n                                  
      do 10581 k=1,n                                                    
      if(abs(u(k)-v(k)) .le. eps)goto 10601                             
      konv=0                                                            
      return                                                            
10601 continue                                                          
10581 continue                                                          
      continue                                                          
      konv=1                                                            
      return                                                            
      end                                                               
      subroutine modstor(iv,av,noth,nz,thv,kpmx,kthmx,istor,fstor,thstor
     *,kp,kth, jerr)
      implicit double precision(a-h,o-z)                                
      double precision thv(nz),fstor(kpmx),thstor(nz,kthmx)             
      integer istor(2,kpmx)                                             
      kp=kp+1                                                           
      jerr=0                                                            
      if(kp .le. kpmx)goto 10621                                        
      jerr=1                                                            
      return                                                            
10621 continue                                                          
      istor(1,kp)=iv                                                    
      fstor(kp)=av                                                      
      if(noth .gt. 0)goto 10641                                         
      istor(2,kp)=0                                                     
      return                                                            
10641 continue                                                          
      kth=kth+1                                                         
      if(kth .le. kthmx)goto 10661                                      
      jerr=2                                                            
      return                                                            
10661 continue                                                          
      istor(2,kp)=kth                                                   
      thstor(:,kth)=thv                                                 
      return                                                            
      end                                                               
      subroutine modpred(no,ni,nz,x,a0,kpl,istor,fstor,thstor,fh)       
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni+nz),fstor(*),thstor(nz,*),fh(no)         
      integer istor(2,*),kpl(2)                                         
      do 10671 i=1,no                                                   
      f=a0                                                              
      do 10681 k=kpl(1),kpl(2)                                          
      is1=istor(1,k)                                                    
      is2=istor(2,k)                                                    
      f=f+fstor(k)*x(i,is1)                                             
      if(is2.eq.0)goto 10681                                            
      xs=x(i,is1)                                                       
      do 10691 l=1,nz                                                   
      f=f+thstor(l,is2)*xs*x(i,l+nz)                                    
10691 continue                                                          
      continue                                                          
10681 continue                                                          
      continue                                                          
      fh(i)=f                                                           
10671 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine modsoln(nz,kpl,istor,fstor,thstor,kv,iv,av,it,kz,tv)   
      implicit double precision(a-h,o-z)                                
      double precision fstor(*), thstor(nz,*),av(*),tv(nz,*)            
      integer kpl(2),istor(2,*),iv(*),it(*)                             
      kv=0                                                              
      kz=0                                                              
      do 10701 k=kpl(1),kpl(2)                                          
      kv=kv+1                                                           
      iv(kv)=istor(1,k)                                                 
      av(kv)=fstor(k)                                                   
      if(istor(2,k) .ne. 0)goto 10721                                   
      it(kv)=0                                                          
      goto 10701                                                        
10721 continue                                                          
      kz=kz+1                                                           
      it(kv)=kz                                                         
      is2=istor(2,k)                                                    
      do 10731 l=1,nz                                                   
      tv(l,kz)=thstor(l,is2)                                            
10731 continue                                                          
      continue                                                          
10701 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
