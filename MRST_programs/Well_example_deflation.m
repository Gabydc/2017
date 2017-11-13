%% Use of Peacemann well models
% In this example we will demonstrate how to set up a flow problems with
% two wells, one rate-controlled, vertical well and one horizontal well
% controlled by bottom-hole pressure. The reservoir is a regular box with
% homogeneous petrophysical properties.
close all
clear all
clc
mrstModule add incomp

%% Set up reservoir model

nf =100;
sz = 20;
szz = 10;
[nxi,nyi,nzi] = deal(1, 1, 1);
[nx,ny,nz] = deal(sz, sz, szz);
[Lx,Ly,Lz] = deal( 500, 500, 25);
rperm = 100;
G = cartGrid([nx,ny,nz],[Lx,Ly,Lz]);
G = computeGeometry(G);
rock = makeRock(G, rperm*milli*darcy, .2);
pl = true;
% Contrast in permeability layers
per = 1;
% Number of layers with same permeability
rlay = 2;
% Create layers of  diverse permeability
if pl == true
    v = [];
    for i = 1:2*rlay:ny
        for j = 0:rlay-1
            if i+j < ny+1
                v = [v j+i];
              
            end
        end
    end


    [I] = Sub2ind_g([1:nx],v,1:nz,nx,ny,nz);
    rock.perm(I) = rperm*10^(-per)*milli*darcy();

end

nf = nf + 1;
f(nf) = figure(nf);
figure(nf)
plotCellData(G, rock.perm/(milli*darcy()),'LineStyle','none'); colorbar
 axis equal tight off


fluid = initSingleFluid('mu', 1*centi*poise,'rho', 1014*kilogram/meter^3);
hT    = computeTrans(G, rock);




%% Add wells wells
W = verticalWell([], G, rock, 1, 1, 1:nz, 'Type', 'bhp', 'Comp_i', 1,...
                'Val', 1.0e8, 'Radius', .12*meter, 'name', 'I');
%disp('Well #1: '); display(W(1));

W = verticalWell(W, G, rock, nx , ny , 1:nz, 'Type', 'bhp', 'Comp_i', 1, ...
            'Val', 1.0e3, 'Radius', .12*meter, 'name', 'P');
%disp('Well #2: '); display(W(2));

state = initState(G, W, 0);



%% Assemble and solve system
gravity reset on;
state1 = incompTPFA(state, G, hT, fluid, 'wells', W, 'MatrixOutput', true);


p0 = state.pressure;
solver = ICCGSolverAD('tolerance', 5.0e-11,'maxIterations', 2000,'cn',0,'x0',p0,'W', W);
%    solver = DPCG_ICSolverAD_cn('tolerance', 5.0e-11,'maxIterations', 2000, 'Z',Z,'cn',0,'x0',p0);
linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
% psolve = @(x,p0,W) ...
%     incompTPFA_Def(x, G, T, fluid, 'wells', W,'LinSolve', linsolve_p);



%state = incompTPFA_Def(state, G, hT, fluid, 'wells', W,'LinSolve', linsolve_p);
A = state1.A;
b = state1.rhs;
[n,m] = size(A); tol = 1e-7;  maxit = 1000; 
if exist('W','var')
    nw = length(W);
    na = n - nw;
    Z = zeros(na,1);
else
     na =n;
     Z = zeros(n,1);
end    

mz = sz/rlay;
for i = 1: mz
   Z(n,i) = 0;
   for j = 1:n/mz
      Z(j+(i-1)*na/mz,i) = 1;
   end
end
size(Z)
Z = sparse(Z); 
L = ichol(A);

solver = DICCGSolverAD('tolerance', 5.0e-11,'maxIterations', 2000,'Z',Z,'cn',0,'x0',p0,'W', W);
%    solver = DPCG_ICSolverAD_cn('tolerance', 5.0e-11,'maxIterations', 2000, 'Z',Z,'cn',0,'x0',p0);
linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
% psolve = @(x,p0,W) ...
%     incompTPFA_Def(x, G, T, fluid, 'wells', W,'LinSolve', linsolve_p);



state = incompTPFA_Def(state, G, hT, fluid, 'wells', W,'LinSolve', linsolve_p);

%   v = diag(A);
%   L = diag(v);
%[L,U] = lu(A);
x0 = zeros(n,1);
%dcg(A,b,Z,tol,maxit,x0)
<<<<<<< HEAD
W = [1];
%[result,flag,res,its,resvec] = ICCG_MRST(A,b,tol,maxit,L,L',x0,'wells',W,'x_true',true);
[result,flag,res,its,resvec] = DICCG_MRST(A,b,Z,tol,maxit,L,L',x0,'wells',W,'x_true',true, 'w_opt', true,'Amatrix_eigs', true);
break
=======
%[result,flag,res,its,resvec] = DICCG_MRST(A,b,Z,tol,maxit,L,L',x0,'wells',W);

%[result,flag,res,its,resvec] = ICCG_MRST(A,b,tol,maxit,L,L',x0,'wells',W);


%%
% We plot the wells to check if the wells are placed as we wanted them.
% (The plot will later be moved to subplot(2,2,1), hence we first find the
% corresponding axes position before generating the handle graphics).
nf = nf + 1;
figure(nf)
subplot(2,2,1), pos = get(gca,'Position'); clf
plotGrid(G, 'FaceColor', 'none');
view(3), camproj perspective, axis tight off,
plotWell(G, W(1), 'radius',  1, 'height', 5, 'color', 'r');
plotWell(G, W(2), 'radius', .5, 'height', 5, 'color', 'b');
>>>>>>> 1915cd1296fc23763d76aa94a9687bd8d0857a31
%% Report results
% We move the plot of the grids and wells to the upper-left subplot. The
% producer inflow profile is shown in the upper-right and the cell
% pressures in the lower-left subplot. In the lower-right subplot, we show
% % the flux intensity, which must be constructed by averaging over cell
% % faces
%figure(nf)
%subplot(2,2,1)
   set(gca, 'Position', pos);  % move the current plot

subplot(2,2,2)
   plot(convertTo(-state.wellSol(2).flux, meter^3/day),'o')
   title('Producer inflow profile [m^3/d]');

subplot(2,2,3)
   plotCellData(G, convertTo(state.pressure(1:G.cells.num), barsa),'EdgeAlpha',.1);
   title('Pressure [bar]')
   view(3), camproj perspective, axis tight off

subplot(2,2,4)
   [i j k] = ind2sub(G.cartDims, 1:G.cells.num);
   I = false(nx,1); I([1 end])=true;
   J = false(ny,1); J(end)=true;
   K = false(nz,1); K([1 end]) = true;
   cf = accumarray(getCellNoFaces(G), ...
      abs(faceFlux2cellFlux(G, state.flux)));
   plotCellData(G, convertTo(cf, meter^3/day), I(i) | J(j) | K(k),'EdgeAlpha',.1);
   title('Flux intensity [m^3/day]')
   view(-40,20), camproj perspective, axis tight, box on
   set(gca,'XTick',[],'YTick',[],'ZTick',[]);

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
