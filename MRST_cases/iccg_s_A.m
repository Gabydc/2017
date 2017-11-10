% This program computes an approximation to x using the deflated 
% preconditioned conjugated gradient method, with splitted preconditioner
function [result,flag,res,its,resvec] = iccg_s_A(A,b,tol,maxit,M1,M2,x0,varargin)
%function diccg_s(n,m,tol,maxit,M1,M2,x0,varargin)
%n = 100; % dimension of matrix
%m = 10;  % number of blocks
% close all
% Select the variables we want to compute
Amatrix_eigs = false;
MAmatrix_eigs = false;
PMAmatrix_eigs = false;

% Select the plots
Plot_residual = true;
Plot_xtrue = true;



%[atype,afun,afcnstr] = iterchk(A);
% if strcmp(atype,'matrix')
%     % Check matrix and right hand side vector inputs have appropriate sizes
     [n,m] = size(A); 
%     if (m ~= n)
%         error(message('MATLAB:dpcg:NonSquareMatrix'));
%     end
%     if ~isequal(size(b),[n,1])
%         error(message('MATLAB:dpcg:RSHsizeMatchCoeffMatrix', m));
%     end
% else
%     n = size(b,1);
%     n = m;
%     if ~iscolumn(b)
%         error(message('MATLAB:dpcg:RSHnotColumn'));
%     end
% end

% Assign default values to unspecified parameters
% if (nargin < 4) || isempty(tol)
%     tol = 1e-6;
% end
% warned = 0;
% if tol < eps
%     warning(message('MATLAB:dpcg:tooSmallTolerance'));
%     warned = 1;
%     tol = eps;
% elseif tol >= 1
%     warning(message('MATLAB:dpcg:tooBigTolerance'));
%     warned = 1;
%     tol = 1-eps;
% end
% if (nargin < 5) || isempty(maxit)
%     maxit = min(n,20);
% end








% Initialize the plots
nf =0;



%[L,U] = lu(A);
%M = eye(n);
i = 1;
x = x0;
%r = P*(b-A*x);
r = b-A*x;
r = M1 \ r;
sr = size(r);
p = M1' \ r;
sp = size(p);
residu(i) = norm(r);
norm(p)
tol = tol* norm(b);
while  (i < maxit) & (residu(i) > tol)
   i = i+1;
   w = A*p;
   alpha = (r'*r)/(p'*w);
   x = x+alpha*p;
   y = M1\w;
   r = r-alpha*y;
   beta = (r'*r)/(residu(i-1)^2);
   z = M1'\r;
   p = z +beta*p;
   residu(i) = norm(r);
  % fout(i) = norm(xtrue-(u+P'*x))/normxtrue;
end

disp(['Number of iterations is: ' num2str(i)])


result = x;
flag = 1;
res = residu(i);
its =i ;
resvec = residu ;






if(Plot_residual == true)
nf = nf+1;
figure(nf)
yplot = residu ;     % plot the relative error
hplot =  semilogy(yplot,'*r');
set(gca,'FontSize', 15)
xlabel('Iteration number')
ylabel('||r_k||_2/||b||_2')
axis tight
end

if(Plot_xtrue == true)
xtrue = A\b;
normxtrue = norm(xtrue);
nf = nf+1;
figure(nf)
vec_x = 1:n;
hplot =  plot(vec_x,x,'*b',vec_x,xtrue,'r');   % plot the solution
set(gca,'FontSize', 15)
xlabel('vec_x')
ylabel('solution')
legend(hplot,'x_k','x_{true}')
axis tight
nf = nf+1;
figure(nf)
hplot =  plot(vec_x,abs(x-xtrue),'*b');   % plot the difference
set(gca,'FontSize', 15)
xlabel('vec_x')
ylabel('abs(x-xtrue)')
axis tight
end
if(Amatrix_eigs == true)
[V,D] = eigs(A,n);
c = cond(A,2);
D = sparse(D);
nf = nf+1;
figure(nf)
% Create axes
axes1 = axes('Parent',figure(nf),'FontSize',14);
for i = 1 : n
plot(V(:,i),'Parent',axes1,'Marker','.','LineStyle','-','color',...
    [0.1*i/(2*n) 0.5*i/(3*n) 0.8*i/n],'DisplayName',['V_' num2str(i)]);
hold on;
end
title(['Eigenvalues A, \kappa (A) = ', num2str(c)],'FontSize',16);
hold off
nf = nf+1;
figure(nf)
% Create axes
axes1 = axes('Parent',figure(nf),'FontSize',14);
for i = 1:n
plot(i,log(D(i,i)),'Parent',axes1,'Marker','*','color',...
    [0.6*i/(2*n) 0.1*i/(6*n) 0.2*i/n],'DisplayName',['\lambda_' num2str(i)]);
hold on;
end
hold off
set(gca,'Xdir','reverse');
xlabel('Eigenvalue','FontSize',16);
ylabel(' log (Value)','FontSize',16);
title('Eigenvalues A','FontSize',16);
end

if(MAmatrix_eigs == true)
    IM = inv(M1);
[V,D] = eigs(IM*A*IM',n);
c = cond(IM*A*IM',2);
D = sparse(D);
nf = nf+1;
figure(nf)
% Create axes
axes1 = axes('Parent',figure(nf),'FontSize',14);
for i = 1 : n
plot(V(:,i),'Parent',axes1,'Marker','.','LineStyle','-','color',...
    [0.1*i/(2*n) 0.5*i/(3*n) 0.8*i/n],'DisplayName',['V_' num2str(i)]);
hold on;
end
title('Eigenvectors M^{-1/2}AM^{-T/2} = L^{-1}AL^{-T}','FontSize',16);
hold off
nf = nf+1;
figure(nf)
% Create axes
axes1 = axes('Parent',figure(nf),'FontSize',14);
for i = 1:n
plot(i,log(D(i,i)),'Parent',axes1,'Marker','*','color',...
    [0.6*i/(2*n) 0.1*i/(6*n) 0.2*i/n],'DisplayName',['\lambda_' num2str(i)]);
hold on;
end
hold off
set(gca,'Xdir','reverse');
xlabel('Eigenvalue','FontSize',16);
ylabel(' log (Value)','FontSize',16);
title({'Eigenvalues M^{-1/2}AM^{-T/2} = L^{-1}AL^{-T}'; ['\kappa (L^{-1}AL^{-T}) = ', num2str(c)]},'FontSize',16);
end




end


