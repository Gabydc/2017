function[p,A,b,h,h1]=defvect_1(p0,nc,cw,q1,div,grad,T,mu,iter,tol,G,nx,ny)

p = initVariablesADI(p0);
 % source term
q = zeros(nc, 1);
 %q(1) = -1.5; q(nc) = 1.5;
 p1=1;
 p2=nx;
 p3=cw;
 p4=nc-nx+1;
 p5=nc; 
   q([p1 p2 p3 p4 p5]) = [q1(1) q1(2) q1(5) q1(4) q1(3)];
%   q(p1)=-1;
%   q(p2)=-1;
%   q(p3)=4;
%   q(p4)=-1;
%   q(p5)=-1;

 eq = - div( (T/mu).*grad(p))+q;
 eq([p1 p2 p3 p4 p5]) = eq([p1 p2 p3 p4 p5])+q([p1 p2 p3 p4 p5]);

l=ichol(eq.jac{1});
[p,h]=ICCG_0(eq.jac{1},eq.val,p0,iter,tol,l);
A=eq.jac{1};
b=eq.val;
figure
h1=plotCellData(G,p);
%figure
% plotCellData(G,A\b);
%  [Va,Da] = eigs(inv(l)*A*inv(l'),length(A));
%  Da=diag(Da);
%  Da=abs(Da);
%  conda=condest(inv(l)*A*inv(l'),2)
%  condadeff=max(Da)/min(Da)
%  figure(40)
% subplot(2,1,1)
% h21=semilogy(Da,'*');
% axis('tight')
% title('M_1A')
%  [Va,Da] = eigs(A,length(A));
%  Da=diag(Da);
%  Da=abs(Da);
%  conda=condest(A,2)
%  condadeff=max(Da)/min(Da)
%   figure(40)
% subplot(2,1,2)
% h21=semilogy(Da,'*');
% title('A')
% axis('tight')
