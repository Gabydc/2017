% This program solves the linear system Ax = b with deflation and split
% preconditioner. The error, true residual and approximation residual are 
% plotted.
%
% programmer: Kees Vuik
% e-mail    : c.vuik@math.tudelft.nl
% date      : 06-10-2005
function dpcg(A,b,Z,tol,maxit,M1,M2,x0,varargin)

%

[n,m] = size(A);
Z  = sparse(Z);
AZ = sparse(A*Z);
E  = Z'*AZ;
EI = inv(E);
P = sparse(eye(n)-AZ*EI*Z');

% deflated cg

xtrue = A\b;
u = Z*EI*Z'*b;
normxtrue = norm(xtrue);

i = 1;
x = x0;
r = P*(b-A*x);
r = M1 \ r;
p = r;
p = M2 \ r;
residu(i) = norm(r);
tresidu(i) = norm(A*x-b);
fout(i)   = norm(xtrue-u)/normxtrue;
tol =tol*norm(b);
while  (i < maxit) && (residu(i) > tol)
   i = i+1;
   PAp = P*(A*p);
   alpha = (r'*r)/(p'*PAp);
   x = x+alpha*p;
   y = PAp;
   y = M1\PAp;
   r = r-alpha*y;
   beta = (r'*r)/(residu(i-1)^2);
   z = p;
   z = M2\p;
   p = r+beta*z;
   residu(i) = norm(r);
   fout(i) = norm(xtrue-(u+P'*x))/normxtrue;
end
xk = (u+P'*x);
%if(Plot_xtrue == true)
xtrue = A\b;
normxtrue = norm(xtrue);
nf = 1;
figure(nf)
vec_x = 1:n;
hplot =  plot(vec_x,xk,'*b',vec_x,xtrue,'r');   % plot the solution
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
%end
% [x]=tdefvec(Z,EI,A,x);
% [qb]=qvec(Z,EI,b);
% x=qb+x;
tr = norm(b-A*(u+P'*x))
disp(['Number of iterations is: ' num2str(i)])
nv = 1 : i; 
yplotre = fout;                  % plot the relative error
yplotrr = residu / norm(u) ;      % plot the relative residual
yplottr = tresidu / norm(u);      % plot the true residual
figure
hplot =  semilogy(nv,yplotre,'*r');
set(gca,'FontSize', 15)
title('Relative Error')
xlabel('Iteration number')
ylabel('(||x-x_k||_2)/||x||_2')
axis tight
figure
%hplot =  semilogy(nv,yplotrr,'*r',nv,yplottr,'*b');
hplot =  semilogy(nv,yplotrr,'*r');
set(gca,'FontSize', 15)
title(['Residual, True Residual ||b-A*x||= ' tr])
%legend('Relative','True')
xlabel('Iteration number')
ylabel('(||r_k||_2)/||b||_2')
axis tight
end


function[qx]=qvec(z,ei,x)
qx=z'*x;
qx=ei*qx;
qx=z*qx;
end
function[px]=defvec(z,ei,a,x)
[qx]=qvec(z,ei,x);
px=x-a*qx;
end
function[px]=tdefvec(z,ei,a,x)
ax=a'*x;
[qax]=qvec(z,ei,ax);
px=x-qax;
end




