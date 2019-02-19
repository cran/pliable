            
      subroutine solveab(nz,b0,grbeta,grtheta,alpha,lam,tt,                     
     1scrat,scrat2,newb,ierrflag)                                         
c solve system of equations for update                                          
c  in pliable lasso procedure                                                   
c it carries out the update in step 4b of robs writeup                          
c model is y=  sum xk*beta(k) +(XkZ)*theta                                      
c  where Z is n by nz  and  XkZ is an n by nz matrix                            
c  formed when each column of Z is multiplied by xk                             
c  this function is called for each predictor k; theta is an nz-vec             
c                                                                               
c inputs                                                                        
c nz: integer- number of modifying variables                                    
c b0: real(nz+1)- initial parameter vector length nz+1                          
c grbeta: real(1)- gradient wrt beta = -sum(xk*r)/n  where xk is X              
c  column of  in current coord descent step, r=current residual                 
c grtheta: real(nz)- gradient wrt theta = -t(XkZ)%*%r/n                         
c alpha: real(1)- mixing par in [0,1]; default .5                               
c lam: real(1)- lasso lam parameter  (on glmnet scale-- ie, refers to           
c  problem  (1/2n)*rss+lam*sum(abs(b))                                          
c  tt: real(1)- backtracking parameter- default 0.1                             
c scrat,scrat2: real(nz) - scratch space                                 
c                                                                               
c  output                                                                       
c     newb: real(nz)- new paramater vector                                      
c  ierrflag: 0 means ok, 1 means solution not found. In calling progrm, exit  loo
c   and set theta vector to 0                                                   
      implicit double precision(a-h,o-z)                                
      integer nz, ierrflag                                                     
      double precision b0(nz+1), grbeta,grtheta(nz),alpha,
     *lam,tt,a(4),b(4)                
      double precision roots(2)                                                           
      double precision val1(4,4),val2(4,4),dng1,dng2,dna,dnb,
     *scrat(nz),scrat2(nz),newb(nz+1)                                           
      
c     write(6,*) 'tt,lam,alpha',tt,lam,alpha                                    
      
      big=10e9                                                                
      eps=1e-3                                                                
      g1=b0(1)-tt*grbeta                                                       
c     write(6,*) 'g1',g1                                                       
c     write(6,*) "b0(2), grtheta(1)",b0(2), grtheta(1)                       
      do i=1,nz                                                                
         scrat(i)=b0(i+1)-tt*grtheta(i)                                        
      end do                                                                 
      
      tt2=tt*alpha*lam                                                         
c     write(6,*) 'tt2',tt2                                                          
      call softthres(scrat,nz,tt2,scrat2)                                      
      
c     write(6,*) 'scrat scrat2',scrat(1),scrat2(1)                            
      
      
      dng1=abs(g1)                                                             
      
      dng2=sqrt(dot_product(scrat2,scrat2))                                    
      
c     write(6,*) 'dng2',dng2                                                         
      
      cc=tt*(1-alpha)*lam                                                     
c     -7.16831982E-02  -2.02317834E-02  -3.59711759E-02         write(6,*) "2*cc   2*
      
      call quadsoln(1d0 ,2*cc, 2*cc*dng2-dng1**2-dng2**2,roots(1),
     1 roots(2))                                                                
c     write(6,*) 'roots',roots(1),roots(2)                                    
      
      a(1)=dng1*roots(1)/(cc+roots(1))                                          
      a(2)=dng1*roots(2)/(cc+roots(2))                                          
      a(3)=dng1*roots(1)/(cc+roots(2))                                          
      a(4)=dng1*roots(2)/(cc+roots(1))                                          
      b(1)=roots(1)*(cc-dng2)/(cc+roots(1))                                     
      b(2)=roots(2)*(cc-dng2)/(cc+roots(2))                                     
      b(3)=roots(1)*(cc-dng2)/(cc+roots(2))                                     
      b(4)=roots(2)*(cc-dng2)/(cc+roots(1))                                     
      
c     write(6,*)'a,b',a,b                                                    
      
      xmin=big
      jhat=1
      khat=1
      do j=1,4                                                                 
         do k=1,4                                                              
            val1(j,k)=big                                                            
            val2(j,k)=big                                                            
            
            den=sqrt(a(j)**2+b(k)**2)                                                 
            
            if(den.gt.0) then                                                        
               val1(j,k)=(1+(cc/den))*a(j)-dng1                                       
               val2(j,k)=(1+cc*(1/b(k) + 1/den))*b(k)-dng2                            
               temp=abs(val1(j,k))+ abs(val2(j,k))                                   
c     write(6,*)  "cc,den,b(k),dng2,xmin,temp",cc,den,b(k),dng2,xmin,temp    
               if(temp.lt.xmin) then                                                 
                  jhat=j                                                                
                  khat=k                                                                
                  xmin=temp                                                             
               end if                                                                
            end if                                                                   
         end do                                                                 
      end do                                                                   
      
      ierrflag=0   
      if(abs(xmin).gt. eps)  then
         ierrflag=1
c     write(6,*) 'jhat,khat,xmin,eps',jhat,khat,xmin,eps 
c     write(6,*) 
c     write(6,*) 'system not solved accurately' 
c     write(6,*) "try rerunning plasso with "
c     write (6,*) "tt= half its current value"
c     write(6,*) "of",tt,";'"
c     write(6,*) "halve the value of tt again"
c     write(6,*) "if this error message persists"
c     write(6,*) ""
         call dblepr('system not solved accurately;',-1,0,0);
         call dblepr('  try rerunning plasso with ',-1,0,0);
         call dblepr(' tt set to half its current value of',-1,tt,1);
         call dblepr(' If error persists,',-1,0,0);
         call dblepr('halve the value of tt again ',-1,0,0);
      end if
      
      dna=a(jhat)                                                              
      dnb=b(khat)                                                              
c     ierrflag=0                                                               
      if((dna.lt.0).or.(dnb.lt.0))  then                                         
c     write(6,*) "one of norm solns negative"  
c     
         ierrflag=1                                                            
      end if                                                                
      xnorm=sqrt(dna**2+dnb**2)                                                  
      
c     write(6,*) 'dna,dnb,xnorm',dna,dnb,xnorm                                          
      
      
c     write(6,*) 'd'                                                                
      
      newb(1)=(b0(1)-tt*grbeta)/(1+cc/xnorm)                                   
c     write(6,*) 'dd'                                                               
      
      do i=1,nz                                                                
         scrat2(i)=b0(i+1)-tt*grtheta(i)                                       
      end do                                                                   
c     write(6,*) 'scrat,tt2',scrat,tt2                                             
      
      
      
      call softthres(scrat2,nz,tt2,newb(2))                                    
      
c     write(6,*) '1 newb',newb                                                      
      
      do i=1,nz                                                                
         newb(i+1)=newb(i+1)/(1+cc*(1/xnorm + 1/dnb))                              
      end do                                                                   
      
c     write(6,*) '2 newb',newb                                                     
      return                                                                   
      end                                                                      
      
      
      
      
      subroutine quadsoln(u,v,w,root1,root2)                                   
c     solve quad eqn u*x**2 + v*x+ w                                               
      implicit double precision(a-h,o-z)      

      double precision u,v,w,root1,root2                                                   
      temp=sqrt(v**2-4*u*w)                                                    
      root1=(-v+temp)/(2*u)                                                    
      root2=(-v-temp)/(2*u)                                                    
      return                                                                   
      end                                                                      
      
      subroutine softthres(x,ni,tt,out)                                        
c     soft threshold components of a  vector                                    
      implicit double precision(a-h,o-z)
      integer  ni                                                               
      double precision x(ni),tt,out(ni)                                                     
      out=0                                                                    
      do i=1,ni                                                                
         if(abs(x(i)).gt.tt) out(i)= sign(1d0,x(i))*(abs(x(i))-tt)              
      end do                                                                   
      return                                                                    
      end                                                                       
