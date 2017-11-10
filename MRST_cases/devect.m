function[p]=defvect(G,N,nc,nf,nx,cw,T,mu)
% G = cartGrid([5 5]);
% Operators
C = sparse([(1:nf )'; (1:nf )'], N, ...
ones(nf,1)*[-1 1], nf, nc);
grad = @(x) C*x;
div = @(x) -C'*x;

% Assemble and solve equations
p = initVariablesADI(zeros(nc,1));
q = zeros(nc, 1);
p1=1;
 p2=nx;
 p3=cw;
 p4=nc-nx+1;
 p5=nc; 
%   q([p1 p2 p3 p4 p5]) = [q1(1) q1(2) q1(5) q1(4) q1(3)];
  q(p1)=-1;
  q(p2)=-1;
  q(p3)=4;
  q(p4)=-1;
  q(p5)=-1;
 % source term
% q(1) = 1; q(nc) = -1;
 % −> quarter five−spot
  eq =  -div( (T/mu).*grad(p))+q;
%eq= div(grad(p))+q;
 % equation
eq(1) = eq(1) + p(1);
 % make solution unique
p = eq.jac{1}\eq.val; % solve equation
plotCellData(G,p);
