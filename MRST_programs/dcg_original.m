% This program solves the linear system Ax = b with deflation and split
% preconditioner. The error, true residual and approximation residual are 
% plotted.
%
% programmer: Gabriela Diaz
% e-mail    : diazcortesgb@gmail.com
% date      : 09-11-2017
function dcg(A,b,Z,tol,maxit,x0,varargin)

Residual = true;
x_true = true;
Convergence = true;
Error = 10^-6;
Amatrix_eigs = true;
MAmatrix_eigs = true;
PMAmatrix_eigs = true;


nf = 0;


[n,m] = size(A);
Z  = sparse(Z);
AZ = sparse(A*Z);
E  = Z'*AZ;
EI = inv(E);
P = sparse(eye(n)-AZ*EI*Z');
%[u]=qvec(Z,EI,b);
u = Z*EI*Z'*b;
if(x_true == true)
    xtrue = A\b;
    normxtrue = norm(xtrue);
end



i = 1;
x = x0;
r = P*(b-A*x);
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
    PAp = P*(A*p);
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
        tresidu(i) = norm(u-A*(u+P'*x));
        fout(i) = norm(xtrue-(u+P'*x))/normxtrue;
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
xk = (u+P'*x);
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
%end
% [x]=tdefvec(Z,EI,A,x);
% [qb]=qvec(Z,EI,b);
% x=qb+x;


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



