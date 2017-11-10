clear all
close all
clc
% Grid and grid information
[nxi,nyi,nzi] = deal(1,1,1);
[nx,ny,nz] = deal(5,5,1);
%[Lx,Ly,Lz] = deal(100,100,1);
%G = cartGrid([nx ny nz],[Lx,Ly,Lz]);
G = cartGrid([nx ny nz]);
G = computeGeometry(G);


%rock.perm = repmat(10*milli*darcy,[G.cells.num, 1]);
rock.perm = repmat(1,[G.cells.num, 1]);
rock.poro = repmat(0.3,[G.cells.num, 1]);


%cr = 1e-6/barsa;
cr=0;
p_r = 200*barsa;
pv_r = poreVolume(G,rock);
pv = @(p)pv_r.*exp(cr*(p-p_r));


mu = 5*centi*poise;
mu=1;
%c = 1e-3/barsa;
c=0;
rho_r = 850*kilogram*meter^3;
rhoS = 750*kilogram/meter^3;
rho = @(p)rho_r.*exp(c*(p-p_r)); 

W = [];
well(1:4)=-1;
    well(5)=4;
wtype    = {'bhp', 'bhp', 'bhp', 'bhp', 'bhp'};
wtarget  = [well(1),   well(2),   well(3),   well(4), well(5)] .* barsa();
wrad     = [0.125, 0.125, 0.125, 0.125, 0.125] .* meter;
wloc     = [  nxi,   nxi,     nx,   nx, ceil(nx/2);
              nyi,   ny,     nyi,   ny, ceil(ny/2)];
wname    = {'W1', 'W2', 'W3', 'W4', 'W5'};
sgn      = [ 1 ,  1 ,  1 ,  1 ,1 ];
     
for w = 1 : numel(wtype),
   W = verticalWell(W, G, rock, wloc(1,w), wloc(2,w), [], ...
                    'Type', wtype{w}, 'Val', wtarget(w), ...
                    'Radius', wrad(w), 'Name', wname{w}, ...
                    'Sign', sgn(w), 'InnerProduct', 'ip_tpf');
end
nc=G.cells.num;
p_init =  zeros(nc, 1);
 q1(1:4)=-0.5;
    q1(5)=2;
    p = initVariablesADI(p_init);
 % source term
q = zeros(nc, 1);
 %q(1) = -1.5; q(nc) = 1.5;
    sz=5;
   cw=(sz*sz)-ceil((sz*sz)/2)+1;
 p1=1;
 p2=nx;
 p3=cw;
 p4=nc-nx+1;
 p5=nc; 
   q([p1 p2 p3 p4 p5]) = [q1(1) q1(2) q1(5) q1(4) q1(3)];


%gravity reset on, g = norm(gravity);

%[z_0,z_max] = deal(0,max(G.cells.centroids(:,3)));
%equil = ode23(@(z,p)g.*rho(p),[z_0,z_max],p_r);
%p_init = reshape(deval(equil,G.cells.centroids(:,3)),[],1);
gravity off



% nc = G.cells.num;
% show=(true(nc,1));
% cellInx=sub2ind(G.cartDims, ...
%     [I-1; I-1; I;  I;  I(1:2)-1], ...
%     [  J;   J; J;  J;  nperf+[2;2]], ... 
%     [K-1;   K; K;  K-1; K(1:2)-[0;1]]);
% show(cellInx) = false;
plotCellData(G,p_init)
plotWell(G,W);


N = G.faces.neighbors;
intInx = all(N ~= 0, 2);
N = N(intInx, :);
hT = computeTrans(G, rock);

cf = G.cells.faces(:,1);
nf = G.faces.num;
T = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]);
T = T(intInx);



 
 
n=size(N,1);
C = sparse([(1:n)'; (1:n)'], N, ...
ones(n,1)*[-1 1], n, nc);
grad = @(x) C*x;
div = @(x) -C'*x;
avg = @(x) 0.5*(x(N(:,1))+x(N(:,2)));

gradz = grad(G.cells.centroids(:,3));

%v = @(p) -(T/mu).*(grad(p)-g*avg(rho(p)).*gradz);
v = @(p) -(T/mu).*(grad(p));
%presEq = @(p) div(avg(rho(p)).*v(p));

wc=W(1).cells;
WI = W(1).WI;
dz = W(1).dZ;
% 
% 
dt =10;
% eq = - div( (T/mu).*grad(p))+q;
presEq = @(p,p0,dt) (1/dt)*(pv(p).*rho(p)-pv(p0).*rho(p0))+div(avg(rho(p)).*v(p));
eq = presEq(p,p_init,dt);
eq([p1 p2 p3 p4 p5]) = eq([p1 p2 p3 p4 p5])+q([p1 p2 p3 p4 p5]);
 A = eq.jac{1};
b = eq.val;
sol = -A\b;

figure
plotCellData(G,sol)
plotWell(G,W);
break
%p_conn =@(bhp) bhp+g*dz.*rho(bhp);
p_conn =@(bhp) bhp;
q_conn = @(p,bhp)WI.*(rho(p(wc))/mu).*(p_conn(bhp)-p(wc));
rateEq = @(p,bhp,qS)qS-sum(q_conn(p,bhp))/rhoS;
ctrlEq = @(bhp)bhp-3;

[p_ad,bhp_ad,qS_ad] = initVariablesADI(p_init,p_init(wc(1)),0);

eq1 = presEq(p_ad);
eq1(wc) = eq1(wc) - q_conn(p_ad,bhp_ad);
eqs = {eq1,rateEq(p_ad,bhp_ad,qS_ad),ctrlEq(bhp_ad)};
eq = cat(eqs{:});
A = eq.jac{1};
b = eq.val;
sol = -A\b;

figure
plotCellData(G,sol(2:length(A)-1))
plotWell(G,W);

break
 eq = - div( (T/mu).*grad(p));
 
p = initVariablesADI(p0);


% Operators





%presEq = @(p,p0,dt)(1/dt)*(pv(p).*rho(p)-pv(p0).*rho(p0)+div(av(rho(p0)).*v(p));


%p_conn =@(bhp) bhp +g*dz.*rho(bhp);



% Assemble and solve equations

q = zeros(nc, 1);
 % source term
q(1) = -1; q(nc) = 1;
 % −> quarter five−spot
eq = -div(grad(p))+q;
 % equation

 % make solution unique
 eq(1) = eq(1) + q(1);
  %eq(nc) = eq(nc) + q(nc);
p = eq.jac{1}\eq.val; % solve equation
plotCellData(G,p);
%l=ichol(eq.jac{1});

A=eq.jac{1};
l=diag(length(A));
b=eq.val;
figure
h1=plotCellData(G,p);
figure
 plotCellData(G,A\b);
 break
 [Va,Da] = eigs(inv(l)*A(2:nx*ny-2,2:nx*ny-2)*inv(l'),length(A)-2);
 Da=diag(Da);
 Da=abs(Da);
 conda=condest(inv(l)*A*inv(l'),2)
 condadeff=max(Da)/min(Da)
 figure(40)
subplot(2,1,1)
h21=semilogy(Da,'*');
axis('tight')
title('M_1A')
 [Va,Da] = eigs(A(2:nx*ny-2,2:nx*ny-2),length(A)-2);
 Da=diag(Da);
 Da=abs(Da);
 conda=condest(A,2)
 condadeff=max(Da)/min(Da)
  figure(40)
subplot(2,1,2)
h21=semilogy(Da,'*');
title('A')
axis('tight')
full(A)