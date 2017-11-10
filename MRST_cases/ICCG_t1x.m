function [xf,ee]=ICCG_t1x(a,b,xi,iteration,tol,l,dir,tolNR)
%a=A;

n=length(a);
x0=xi;
r0=b-a*xi;
lb=l\b;
r0=l\r0;
p0=l'\r0;
nor=abs(lb'*lb);

%nor=1;    
text = [dir 'error.txt'];
fileID = fopen(text,'a');
fprintf(fileID,'Iteration & $||M^{-1}b||$ & $||M^{-1}r||$ & $\frac{||M^{-1}r||}{||M^{-1}b||}$ & $max(abs(b-Ax))$ \\\\ \n');
fclose(fileID);
ax=norm(x0);
if norm(x0) <= tolNR
    xf=x0;
    ee=norm(x0);
    return
end
for iter=1:iteration 

    ap=a*p0;
    alpha=(r0'*r0)/((ap)'*p0);
    xf=xi+alpha*p0;
    r=r0-alpha*(l\ap);
    beta=(r'*r)/(r0'*r0);
    p=l'\r+beta*p0;
    p0=p;
    r0=r;
    rr=abs(r'*r);
ee(iter)=abs(r'*r)/nor;
erb=ee(iter);
 e=max(abs(b-a*xf));
fileID = fopen(text,'a');

fprintf(fileID,'%d & %1.1e & %1.1e & %1.1e & %1.1e \\\\ \n', iter,nor,rr,erb,e);
fclose(fileID);  
    %pause
    flag=0;
    xif=(xi-xf);
    ex=abs(xif'*xif)/abs(xi'*xi);
    if (ee>=tol)
        flag=1;
    end
    if flag==0
        break
    end
    xi=xf;

end
end