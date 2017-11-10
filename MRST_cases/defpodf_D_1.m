%This program gives the matrix for the solution of the Laplace problem, the
%linear system is solved with the CGF conjugate gradient function

function [U,S,hd]=defpodf_D(x,p)
lim=0.1*p;
%x=zd;
lx=size(x,1);
ly=size(x,2);
D=x'*x;
[V,L] = eigs(D,ly);
S=zeros(lx,ly);
g=diag(L);
g1=sqrt(g);
g1=1./g1;

%[V,D] = eigs(X,n);
figure(3000)
hd=plot(g1,'og'); 
%axis('tight')
title('Eigenvalues X');
ylabel('Value ')
xlabel('Eigenvalue')
gmax=max(g1);
s=0;
for i=1:ly
    if g1(i)>lim*gmax
        s=s+1;
        gm(s)=g1(i);
    end
end
S(1:s,1:s)=diag(gm);
U=x*V*S';
figure(3001)
hd=plot(gm,'og'); 
%axis('tight')
title('Eigenvalues X');
ylabel('Value ')
xlabel('Eigenvalue')
% file='eig_pod';
% B=[dir '/results/qfs_1/' file '.fig'];
% saveas(h1,B)
% B=[dir '/results/qfs_1/' file '.jpg'];
% saveas(h1,B)
