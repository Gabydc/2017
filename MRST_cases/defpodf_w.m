%This program gives the matrix for the solution of the Laplace problem, the
%linear system is solved with the CGF conjugate gradient function

function [V,D,hd]=defpodf_w(x,n)

X=x*x';
[V,D] = eigs(X,n);
figure(3000)
hd=plot(diag(D),'*r'); 
axis('tight')
hold on
title('Eigenvalues X');
ylabel('Value ')
xlabel('Eigenvalue')


% file='eig_pod';
% B=[dir '/results/qfs_1/' file '.fig'];
% saveas(h1,B)
% B=[dir '/results/qfs_1/' file '.jpg'];
% saveas(h1,B)
