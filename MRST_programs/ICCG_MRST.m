% This program solves the linear system Ax = b with deflation and split
% preconditioner. The error, true residual and approximation residual are
% plotted.
%
% programmer: Gabriela Diaz
% e-mail    : diazcortesgb@gmail.com
% date      : 09-11-2017
function [result,flag,res,its,resvec] = ICCG_MRST(A,b,tol,maxit,M1,M2,x0,varargin)
warning on backtrace
%warning off verbose


opt = struct( 'Residual', false, ...
'x_true', false, ...
'Convergence' , false, ...
'Amatrix_eigs', false, ...
'MAmatrix_eigs', false, ...
'Iter_m' , false, ...
'A_cn' , false, ...
'dir', [], ...
'Error', '10^-6', ...
'nf', '1000', ...
'wells', []);
 opt = merge_options(opt, varargin{:});
 Error = opt.Error;
 nf = opt.nf;

 dir = opt.dir;
 W  = opt.wells;
display(W)
[n,m] = size(A);
nw = numel(W);
if nw > 0
    na = n - nw;
else
     na = n;
end 
if (opt.A_cn )
    SZ_A =size(A);
    C_A = cond(A);
    disp(['The size of A is: ' num2str(SZ_A)])
    disp(['The condition number of A is: ' num2str(C_A)])
end

A(na+1 : n,:);
b(na+1 : n);
A(:,na+1 : n);
if (m ~= n)
    warning(['NonSquareMatrix'])
    error(message());
end
if ~isequal(size(b),[n,1])
   warning('RSHsizeMatchCoeffMatrix', m);
end
if (nargin < 4) || isempty(tol)
    tol = 1e-6;
end
warned = 0;
if tol < eps
    warning('tooSmallTolerance');
    warned = 1;
    tol = eps;
elseif tol >= 1
    warning(['tooBigTolerance']);
    warned = 1;
    tol = 1-eps;
end
if (nargin < 5) || isempty(maxit)
    maxit = min(n,20);
end

if (nw > 0)
for i = 1 : nw
    i;
    n-nw+i;
    bw = b(n-nw+i,1);
    Aw = A(n-nw+i,n-nw+i);
pw(n-nw+i,1) = bw / Aw;
end

A1(1:n-nw,1:n-nw)=A(1:n-nw,1:n-nw);
M11(1:n-nw,1:n-nw)=M1(1:n-nw,1:n-nw);
M22(1:n-nw,1:n-nw)=M2(1:n-nw,1:n-nw);
b1(1:n-nw,1)=b(1:n-nw,1);
x01(1:n-nw,1)=x0(1:n-nw,1);
A = A1;
M1 = M11;
M2 = M22;
b = b1;
x0 = x01;
clear A1 Z1 M11 M22 b1 x01
end




size(b)
size(A)
size(x0)



if(opt.x_true )
    xtrue = A\b;
    normxtrue = norm(xtrue);
end



i = 1;
x = x0;
Mb = M1 \ b;
r = b - A * x;
r = M1 \ r;
%p = r;
p = M2 \ r;
residu(i) = norm(r);
if(opt.x_true )
    fout(i)   = norm(xtrue-x)/normxtrue;
    tresidu(i) = norm(b-A*x);
end

tol =tol*norm(Mb);
if(residu(i) < tol) 
   disp(['DICCG only needs one iteration, initial residual is, r_1 = ||P*M^{-1}(b-A*xk)||_2' num2str(residu(i))])
   disp(['True residual is, r_1 = ||b-A*xk||_2: ' num2str(tresidu(i))])
end
while  (i < maxit) && (residu(i) > tol)
    
    if(opt.Convergence) && (residu(i) < Error)
        % If the residual increases, the approximation will be the previous
        % solution
        xacc = x;
    end
    
    i = i+1;
    w = A * p;
    alpha = (r' * r) / (p' * w);
    x = x + alpha * p;
    y = M1 \ w;
    r = r - alpha * y; 
    beta = (r' * r)/(residu(i-1)^2);
    %z = r;
     z = M2 \ r;
    p = z + beta * p;
    residu(i) = norm(r);
    if(opt.x_true )
        xk = x;
        tresidu(i) = norm(b-A*xk);
        fout(i) = norm(xtrue-xk)/normxtrue;
    end
    
    if(opt.Convergence ) && (residu(i) < Error)
        % If the residual increases, the approximation will be the previous
        % solution
        if (residu(i) >= residu(i-1))
            flag = 0;
        end
        if flag == 0
            warning(['Maximum accuracy is : ' num2str(residu(i))])
            break
        end
    end
    
end
if (opt.Iter_m )
disp(['Number of iterations is: ' num2str(i)])
end

xk = x;
tr = norm(b-A*xk)/norm(b);
Mr = M1 \ (b-A*xk);
ptr = norm(Mr)/norm(Mb);


result = xk;
flag= 0;
res = residu(i);
its =i;
resvec = residu;

if (nw > 0)
    for i = 1 : nw
        result(n-nw+i,1) = pw(n-nw+i,1);
    end
end

if(opt.x_true )
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
if(opt.x_true )
    yplotre = fout;                  % plot the relative error
    yplottr = tresidu / norm(Mb);      % plot the true residual
    nf = nf+1;
    figure(nf)
    hplot =  semilogy(nv,yplotre,'*r');
    set(gca,'FontSize', 16)
    title('Relative Error')
    xlabel('Iteration number','FontSize', 16)
    ylabel('(||x-x_k||_2)/||x||_2','FontSize', 16)
    axis tight
end
if(opt.Residual )
    yplotrr = residu / norm(Mb) ;      % plot the relative residual
    nf = nf+1;
    figure(nf)
    %hplot =  semilogy(nv,yplotrr,'*r',nv,yplottr,'*b');
    hplot =  semilogy(nv,yplotrr,'*r');
    set(gca,'FontSize', 12)
    title({['Relative Residual']; ['||b-A*x_k||_2/||b||_2= ' num2str(tr)];...
        ['||M^{-1}(b-A*x_k)||_2/||M^{-1}b||_2= ' num2str(ptr)]})
    %legend('Relative','True')
    xlabel('Iteration number','FontSize', 16)
    ylabel('||M^{-1}r_k||_2/||M^{-1}b||_2','FontSize', 16)
    axis tight
end

 
%% Eigenvalues eigenvectors
if(opt.Amatrix_eigs )
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

if(opt.MAmatrix_eigs )
    
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





end

