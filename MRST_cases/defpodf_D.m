%This program gives the matrix for the solution of the Laplace problem, the
%linear system is solved with the CGF conjugate gradient function

function [U,S,hd]=defpodf_D(x,dir)
%x=zd;
lx=size(x,1);
ly=size(x,2);
D=x'*x;
[V,L] = eigs(D,ly);
S=zeros(lx,ly);
g=diag(L);
g1=sqrt(g.*g);
g1=1./g1;
S(1:ly,1:ly)=diag(g1);
S=sparse(S);
V=sparse(V);
x=sparse(x);
U=x*V*S';
%[V,D] = eigs(X,n);
figure(3000)
hd=plot(log(g1),'ob'); 
axis('tight')
%title('Eigenvalues X');
ylabel('log(Value) ')
xlabel('Eigenvalue')
hold on


file='eig_pod';
B=[dir file '.fig'];
saveas(hd,B)
B=[dir  file '.jpg'];
saveas(hd,B)
