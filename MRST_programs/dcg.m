% This program solves the linear system Ax = b with deflation and split
% preconditioner. The error, true residual and approximation residual are
% plotted.
%
% programmer: Gabriela Diaz
% e-mail    : diazcortesgb@gmail.com
% date      : 09-11-2017
function dcg(A,b,Z,tol,maxit,x0,varargin)

Residual = false;
x_true = false;
Convergence =false;
Error = 10^-6;
Amatrix_eigs = true;
MAmatrix_eigs = true;
PMAmatrix_eigs = true;
nf = 0;


[n,m] = size(A);
if (m ~= n)
    error(message('MATLAB:dcg:NonSquareMatrix'));
end
if ~isequal(size(b),[n,1])
    error(message('MATLAB:dcg:RSHsizeMatchCoeffMatrix', m));
end
[nz,mz] = size(Z);
if ~isequal(size(b),[nz,1])
    error(message('MATLAB:dcg:WrongDeflationMAtrixSize', nz));
end

if (nargin < 4) || isempty(tol)
    tol = 1e-6;
end
warned = 0;
if tol < eps
    warning(message('MATLAB:dcg:tooSmallTolerance'));
    warned = 1;
    tol = eps;
elseif tol >= 1
    warning(message('MATLAB:dcg:tooBigTolerance'));
    warned = 1;
    tol = 1-eps;
end
if (nargin < 5) || isempty(maxit)
    maxit = min(n,20);
end








Z  = sparse(Z);
AZ = sparse(A*Z);
E  = Z'*AZ;
EI = inv(E);
%P = sparse(eye(n)-AZ*EI*Z');
[u]=qvec(Z,EI,b);
%u = Z*EI*Z'*b;
if(x_true == true)
    xtrue = A\b;
    normxtrue = norm(xtrue);
end



i = 1;
x = x0;
r = b - A * x;
r = defvec(Z,EI,A,r);
%r = P * r;
%r = M1 \ r;
p = r;
%p = M2 \ r;
residu(i) = norm(r);
if(x_true == true)
    fout(i)   = norm(xtrue-u)/normxtrue;
end

tol =tol*norm(b);
while  (i < maxit) && (residu(i) > tol)
    
    if(Convergence == true) && (residu(i) < Error)
        % If the residual increases, the approximation will be the previous
        % solution
        xacc = x;
    end
    
    i = i+1;
    w = A * p;
    %PAp = P*(A*p);
    PAp = defvec(Z,EI,A,w);
    alpha = (r'*r)/(p'*PAp);
    x = x+alpha*p;
    y = PAp;
    %  y = M1\PAp;
    r = r-alpha*y;
    beta = (r'*r)/(residu(i-1)^2);
    z = p;
    % z = M2\p;
    p = r+beta*z;
    residu(i) = norm(r);
    if(x_true == true)
        % tresidu(i) = norm(u-A*(u+P'*x));
        [xk]=tdefvec(Z,EI,A,x);
        [Qb] = qvec(Z,EI,b);
        xk = Qb + xk;
        tresidu(i) = norm(b-A*xk);
        % fout(i) = norm(xtrue-(u+P'*x))/normxtrue;
        fout(i) = norm(xtrue-xk)/normxtrue;
    end
    
    if(Convergence == true) && (residu(i) < Error)
        % If the residual increases, the approximation will be the previous
        % solution
        if (residu(i) >= residu(i-1))
            flag = 0;
        end
        if flag == 0
            disp(['Maximum accuracy is : ' num2str(residu(i))])
            break
        end
    end
    
end
%xk = (u+P'*x);
[xk] = tdefvec(Z,EI,A,x);
[Qb] = qvec(Z,EI,b);
xk = Qb + xk;
tr = norm(b-A*xk)/norm(b);

disp(['Number of iterations is: ' num2str(i)])
if(x_true == true)
    normxtrue = norm(xtrue);
    nf = nf + 1;
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
    hplot =  plot(vec_x,abs(xk-xtrue),'*b');   % plot the difference
    set(gca,'FontSize', 15)
    xlabel('vec_x')
    ylabel('abs(x-xtrue)')
    axis tight
end

[n,m] = size(A);
nv = 1 : i;
if(x_true == true)
    yplotre = fout;                  % plot the relative error
    yplottr = tresidu / norm(b);      % plot the true residual
    nf = nf+1;
    figure(nf)
    hplot =  semilogy(nv,yplotre,'*r');
    set(gca,'FontSize', 15)
    title('Relative Error')
    xlabel('Iteration number')
    ylabel('(||x-x_k||_2)/||x||_2')
    axis tight
end
if(Residual == true)
    yplotrr = residu / norm(b) ;      % plot the relative residual
    nf = nf+1;
    figure(nf)
    %hplot =  semilogy(nv,yplotrr,'*r',nv,yplottr,'*b');
    hplot =  semilogy(nv,yplotrr,'*r');
    set(gca,'FontSize', 15)
    title({['Residual, True Residual']; ['||b-A*x_k||= ' num2str(tr)]})
    %legend('Relative','True')
    xlabel('Iteration number')
    ylabel('(||b-A*x_k||_2)/||b||_2')
    axis tight
end

 
%% Eigenvalues eigenvectors
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
title('Eigenvectors A','FontSize',16);
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
title(['Eigenvalues A, \kappa (A) = ', num2str(c)],'FontSize',16);
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
 [ymax] = max(diag(D));
 [ymin] = min(diag(D));
 ylim(axes1,[log(ymin) log(ymax)])
set(gca,'Xdir','reverse');
xlabel('Eigenvalue','FontSize',16);
ylabel(' log (Value)','FontSize',16);
title({'Eigenvalues M^{-1/2}AM^{-T/2} = L^{-1}AL^{-T}'; ['\kappa (L^{-1}AL^{-T}) = ', num2str(c)]},'FontSize',16);

end

if(PMAmatrix_eigs == true)
    
    
     IM = inv(M1);
    Q = Z * EI * Z';
    P = sparse(eye(n)-AZ*EI*Z');
[V,D] = eigs(P*IM*A*IM',n);
D = diag(D);
D = real(D(mz+1:n));
lmax = max(D);
lmin = min(D);
c = lmax/lmin;

nf = nf+1;
figure(nf)
% Create axes
axes1 = axes('Parent',figure(nf),'FontSize',14);
for i = mz + 1 : n
plot(V(:,i),'Parent',axes1,'Marker','.','LineStyle','-','color',...
    [0.1*i/(2*n) 0.5*i/(3*n) 0.8*i/n],'DisplayName',['V_' num2str(i)]);
hold on;
end
title('Eigenvectors PM^{-1/2}AM^{-T/2} = PL^{-1}AL^{-T}','FontSize',16);
hold off
nf = nf+1;
figure(nf)
% Create axes
axes1 = axes('Parent',figure(nf),'FontSize',14);
for i = 1 : n -mz
plot(i,log(D(i)),'Parent',axes1,'Marker','*','color',...
    [0.6*i/(2*n) 0.1*i/(6*n) 0.2*i/n],'DisplayName',['\lambda_' num2str(i)]);
hold on;
end
hold off
ylim(axes1,[log(ymin) log(ymax)])
set(gca,'Xdir','reverse');
xlabel('Eigenvalue','FontSize',16);
ylabel(' log (Value)','FontSize',16);
title({'Eigenvalues PM^{-1/2}AM^{-T/2} = PL^{-1}AL^{-T}'; ['\kappa_{eff} (PL^{-1}AL^{-T}) = ', num2str(c)]},'FontSize',16);
end

if(PMAmatrix_eigs == true)
    
     IM = inv(M1);
    Q = Z * EI * Z';
    P = sparse(eye(n)-AZ*EI*Z');
[V,D] = eigs(P*IM*A*IM',n);
D = diag(D);
D = real(D);
lmax = max(D);
lmin = min(D);
c = lmax/lmin;
nf = nf+1;
figure(nf)
% Create axes
axes1 = axes('Parent',figure(nf),'FontSize',16);
for i =  1 : n
plot(V(:,i),'Parent',axes1,'Marker','.','LineStyle','-','color',...
    [0.1*i/(2*n) 0.5*i/(3*n) 0.8*i/n],'DisplayName',['V_' num2str(i)]);
hold on;
end
title('Eigenvectors PM^{-1/2}AM^{-T/2} = PL^{-1}AL^{-T}, all','FontSize',16);
hold off
nf = nf+1;
figure(nf)
% Create axes
axes1 = axes('Parent',figure(nf),'FontSize',16);
for i =  1 : n
plot(i,log(D(i)),'Parent',axes1,'Marker','*','color',...
    [0.6*i/(2*n) 0.1*i/(6*n) 0.2*i/n],'DisplayName',['\lambda_' num2str(i)]);
hold on;
end
hold off
set(gca,'Xdir','reverse');
xlabel('Eigenvalue','FontSize',16);
ylabel(' log (Value)','FontSize',16);
title({'Eigenvalues PM^{-1/2}AM^{-T/2} = PL^{-1}AL^{-T}, all'; ['\kappa_{eff} (PL^{-1}AL^{-T}) = ', num2str(c)]},'FontSize',16);
end
















end



function[Qx]=qvec(Z,EI,x)
Qx=Z'*x;
Qx=EI*Qx;
Qx=Z*Qx;
end
function[Px]=defvec(Z,EI,A,x)
[Qx]=qvec(Z,EI,x);
Px=x-A*Qx;
end
function[Px]=tdefvec(Z,EI,A,x)
Ax=A'*x;
[QAx]=qvec(Z,EI,Ax);
Px=x-QAx;
end