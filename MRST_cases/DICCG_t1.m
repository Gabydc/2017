
function[xf,ee]=DICCG_t1(a,b,xi,iteration,tol,z,l)

n=size(a,2);
e=z'*a*z;
ei=sparse(inv(e));
%size(ez)


[pb]=deflatevect(z,ei,a,b);
lb=l\b;
plb=l\pb;
nor1=abs(lb'*lb);
r0=b-a*xi;
[r0]=deflatevect(z,ei,a,r0);
r0=l\r0;
p0=l'\r0;
nor=abs(lb'*lb);
%nor=1;
for iter=1:iteration
         [ap]=deflatevect(z,ei,a,a*p0);
         [apt]=tdeflatevect(z,ei,a,p0);
     alpha=(r0'*r0)/((ap)'*p0);   
     xf=xi+alpha*p0;
     r=r0-alpha*(l\ap);
     beta=(r'*r)/(r0'*r0);
     p=l'\r+beta*p0;  
     p0=p;
     r0=r;
     color=[0.3 0.5 0.7];

      ee(iter)=abs(r'*r)/nor;

       flag=0;

     if (ee(iter)>=tol)
         flag=1;
     end
     if flag==0
         break
     end     
     xi=xf;
end
[xf]=tdeflatevect(z,ei,a,xf);
qb=z'*b;
qb=ei*qb;
qb=z*qb;
xf=qb+xf;