function [xf,flag,relres,iter,resvec,conda] = pcg_g(A,b,tol,maxit,M1,M2,x0,varargin)
size(b)
size(A)

x0 = zeros(size(b));
n=size(A,2);
r0=b-A*x0

r0=M1\r0;
p0=M2\r0;
conda=condest(inv(M1*M2)*A);
[vach,dach]=eigs(inv(M1*M2)*A);

for iter=1:maxit
     w = A*p0;
     rn0 = r0'*r0;
     alpha = rn0/(w'*p0);
     xf=x0+alpha*p0
     r=r0-alpha*(M1\w);
     rn = r'*r; 
     beta=rn/rn0;
     p=(M2\r)+beta*p0 ;
     p0=p;
     r0=r;
    if norm(xf)==0
       e=0;
       resvec(iter) = norm(r);
     else
          e=abs((xf-x0)./xf)*100;
          resvec(iter) = norm(r);
     ee=sqrt(e'*A*e);
     
     color=[0.2 0.8 0.6];
     figure(123)
     hline=plot(iter,log(ee),'o','Color',color);
     hold on
     hline=plot(iter,norm(r),'*','Color','b');
    end
    x0=xf;
       flag=0;
     if (ee>=tol)
         flag=1;
     end
     if flag==0
         break
     end     
     
    
     
 end