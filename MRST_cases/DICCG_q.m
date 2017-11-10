
function[xf,iter,ee,hl1,h2,h21]=DICCG_q(a,b,xi,iteration,tol,z,l,eig,text,nn,dir,def,fc)
hl=0;
hl1=0;
hl11=0;
%tol=10^-5;
z=sparse(z);
n=size(a,2);
e=z'*a*z;
ie=inv(e);

if eig==1
n=size(a,2);
fprintf('DICCG')
[pl]=deflatevect(z,ie,a,a*inv(l'));
[Va,Da] = eigs(inv(l)*pl,n);
 Da=diag(Da);
 Da=abs(Da);
 conda=condest(inv(l)*pl,2)
 condadeff=max(Da)/min(Da)
 
 fileID = fopen(text,'a');
fprintf(fileID,'%6s %6s %6.2e %6s\n','DCGCh',' &',condadeff, ' \\');
fclose(fileID);
 text1 = [dir def 'eigs_DCGCh.txt'];
 fileID = fopen(text1,'w');
fprintf(fileID,'%6.2e \n', Da);
fclose(fileID);
file=[dir def 'eigs_DCGCh'];
 save(file,'Da')
figure(3)
subplot(2,1,2)
h2=semilogy(Da,'*');
axis('tight')
title('DICCG')
ylabel('log scale)')
xlabel('Eigenvalue')
hold on
figure(4)
subplot(2,1,2)
h21=semilogy(Da,'*');
axis([1, n, 0.004, 1.5])
title('DICCG')
ylabel('log scale')
xlabel('Eigenvalue')
hold on
else h2=0; h21=0;
end

%size(ez)
ez=ie*z';
q=z*ez;
pd=eye(n)-a*q;
  pda=pd*a;
  pb=pd*b;
% pb=z'*b;
% pb=ie*pb;
% pb=z*pb;
% pb=a*pb;
% pb=b-pb;
%[pb]=deflatevect(z,ie,a,b);
lb=l\b;
plb=l\pb;
nor1=abs(plb'*plb);
r0=b-a*xi;
%[r0]=deflatevect(z,ie,a,r0);
r0=pd*r0;
r0=l\r0;
p0=l'\r0;
nor=abs(lb'*lb);
%nor=1;
for iter=1:iteration
% ap=a*p0;
% ap=z'*ap;
% ap=ie*ap;
% ap=z*ap;
% ap=a*ap;
% ap=p0-ap;
%[ap]=deflatevect(z,ie,a,p0);
         ap=pda*p0; 
     alpha=(r0'*r0)/((ap)'*p0);   
     xf=xi+alpha*p0;
     r=r0-alpha*(l\ap);
     beta=(r'*r)/(r0'*r0);
     p=l'\r+beta*p0;  
     p0=p;
     r0=r;
     color=[0.1 0.5 0.5];
%     if xf==0
%        e=0;
%      else
%           e=(abs(xi-xf)./abs(xf))*100;
%      ee=sqrt(e'*pda*e);
%      xtrue=q*b+pd'*xi;
%      rt=(b-a*xtrue);
%      rt=abs(rt'*rt)/nor1;
%       figure(1)
%       hl=semilogy(iter,(rt),'p','Color',color);
%       hold on
% %     end
%     
%       fra=abs(r);
%      b1=abs(b);
%       for rr=1:n
%           if b1(rr)==0
%               e1(rr)=0;
%           else             
%               e1(rr)=fra(rr)/b1(rr);
%           end
%       end
%      ee1=sqrt(e1*pda*e1');
%      xff=q*b+pd'*xf;
%      xfff=(a*xff-b)'*(a*xff-b);
      ee=abs(r'*r)/nor;
     figure(fc)
     hl1=semilogy(iter,ee,'p','Color',color);
     hold on
     
%       if xf==0
%          e11=0;
%      else
%           e11=(sqrt((xi-xf)'*pda*(xi-xf))./sqrt((xf')*pda*(xf)))*100;
%      ee11=sqrt(e'*pda*e);
% ed=pb-pda*xi;
% ed=abs(ed'*ed)/nor1;
%          figure(6)
% 
%       hl11=semilogy(iter,(ed),'p','Color',color);
%      
%       hold on
% %      end
     
       flag=0;

     if (ee>=tol)
         flag=1;
     end
     if flag==0
         break
     end     
     xi=xf;
end
% xpt=a'*xf;
% xpt=z'*xpt;
% xpt=ie*xpt;
% xpt=z*xpt;
% xpt=xf-xpt;
% [xpt]=tdeflatevect(z,ie,a,xf);
% qb=z'*b;
% qb=ie*qb;
% qb=z*qb;
%xf=qb+xpt;
xf=q*b+pd*xf;