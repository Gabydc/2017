%% Single phase flow simulation using AD
% This example goes through the steps of setting up a single-phase
% simulation with a single horizontal well using the automatic
% differentiation framework.

clf
% Required modules
try
   require ad-fi
catch
   mrstModule add ad-fi
end

% Setup 10x10x10 grid of 200x200x50 m model.
nx = 35;    ny = 35;    nz = 1;
Dx = 200;   Dy = 200;   Dz = 50;
G = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
%G = cartGrid([nx, ny, nz]);
G = computeGeometry(G);

%% Setup rock properties
%

% Assume homogeneous/isotropic rock.
permX = 1*milli*darcy;
poro  = 0.3;
rock.perm = repmat(permX, [G.cells.num, 1]);
rock.poro = repmat(poro , [G.cells.num, 1]);
 sz=nx;
 lsz=round(sz/5); 
per=1;
for i= 1:2:5
    rock.perm(1+(i-1)*lsz*nx:i*lsz*nx)  = 10^(-per)*milli*darcy();
end
% for i=2:2:6
%  %for i=2
%  if i==2
%      init1=lsz*(i-i/4);
%   rock.perm(1+init1:init1+lsz)  = repmat(10^(-per)*milli*darcy(), [lsz, 1]);
%  else
%      init=init1+lsz*(i-2);
%      rock.perm(init+1:init+lsz)  = repmat(10^(-per)*milli*darcy(), [lsz, 1]);
%  end
%  end
% for i=1:2:8
%  rock.perm(1+lsz*(i-1):lsz*i)  = repmat(10^(-per)*milli*darcy(), [lsz, 1]);
% end

% Set rock compressibility:
cr = 1e-6/barsa;

%%
% In the case of non-zero rock-compressiblity cr, the input rock porosity
% is taken as reference at a given reference pressure p_r. The grid pore
% volume (pv) becomes a function of pressure given by the differential
% equation
%                 cr = (d pv/d p)/pv
% which results in
%                 pv(p) = pv_r e^( cr(p-p_r) )
% in which pv_r is the reference pore volume (rock.poro x volume) and p_r
% is the reference pressure. We assume the reference pressure is 200 Bar:

pv_r = poreVolume(G, rock);
p_r  = 200*barsa;

% Finally, the pressure-dependent function for pore-volumes becomes:

pv   = @(p) pv_r .* exp( cr .* (p - p_r) );

%% Fluid (oil) properties
% We Assume constant viscosity:
mu   = 5*centi*poise;
% Assume that the oil compressibility can be approximated as constant in
% the reservoir:
c    = 1e-3/barsa;
% When the compressibility is constant, the fluid density becomes a
% function of pressure given by the differential equation
%                 c = (d rho/d p)/rho
% which results in
%                 rho(p) = rho_r e^( c(p-p_r) )
% in which rho_r is the reference fluid density at reference pressure p_r.
% We assume that rho_r = 850 kg/m^3 at p_r = 200 Bar:

rho_r = 850*kilogram/meter^3;

% Finally define the pressure dependent function for rho:
rho   = @(p) rho_r .* exp( c .* (p - p_r) );

% Computing surface volume rates requires fluid density at surface
% conditions. We assume it is 750 kg/m^3:
rhoS = 750*kilogram/meter^3;

%% Single horizontal well in J direction
% We consider a well with 8 connections in the J direction. 

W = [];
nperf = 1;
I = repmat(1, [nperf, 1]);
J = 1;
K = repmat(1, [nperf, 1]);

% Convert IJK-indices to linear index (as used in G)
cellInx = sub2ind(G.cartDims, I, J, K);

W = addWell(W, G, rock, cellInx, 'Name', 'producer', 'InnerProduct', 'ip_tpf', 'sign', -1);

nperf = 1;
I = repmat(nx, [nperf, 1]);
J = 1; %(1 : nperf).';
K = repmat(1, [nperf, 1]);

cellInx = sub2ind(G.cartDims, I, J, K);

W = addWell(W, G, rock, cellInx, 'Name', 'producer', 'InnerProduct', 'ip_tpf', 'sign', -1);

nperf = 1;
I = 1;
J = ny; %(1 : nperf).';
K = repmat(1, [nperf, 1]);

% Convert IJK-indices to linear index (as used in G)
cellInx2 = sub2ind(G.cartDims, I, J, K);

W = addWell(W, G, rock, cellInx2, 'Name', 'producer', 'InnerProduct', 'ip_tpf', 'sign', 1);

nperf = 1;
I = nx;
J = ny;
K = repmat(1, [nperf, 1]);

% Convert IJK-indices to linear index (as used in G)
cellInx = sub2ind(G.cartDims, I, J, K);

W = addWell(W, G, rock, cellInx, 'Name', 'producer', 'InnerProduct', 'ip_tpf', 'sign', -1);

nperf = 1;
I = ceil(nx/2);
J = ceil(ny/2);
K = repmat(1, [nperf, 1]);

% Convert IJK-indices to linear index (as used in G)
cellInx = sub2ind(G.cartDims, I, J, K);

W = addWell(W, G, rock, cellInx, 'Name', 'injector', 'InnerProduct', 'ip_tpf', 'sign', -1);

% Plotting
f = [figure(1), figure(2)];
set(0, 'CurrentFigure', f(1));
clf
 hp=plotCellData(G, rock.perm);
%  
%  view(0,90) 
%plotGrid(G, 'FaceColor', 'g', 'FaceAlpha', .3, 'EdgeColor', 'w');
plotWell(G, W);
axis off;
set(f(1), 'Color', 'w'); 
camproj perspective;
view(3);
%break

%% Initial conditions
% We assume that the reservoir is initially at equilibrium. This means that
% the following condition must be satisfied:
%          dp/dz = g rho,
% where g is the gavitational accelleration. This relation can be solved
% analytically for p, but alternatively one can solve the above ODE with
% 'initial condtition' p(z_0) = p_r:

gravity reset on
g = norm(gravity);
[z_0, z_max] = deal(0, max(G.cells.centroids(:,3)) + 5);
press  = ode23(@(z,p) g .* rho(p), [z_0, z_max], p_r);
% We then interpolate onto the grid using cell centers:
p_init = reshape(deval(press, G.cells.centroids(:,3)), [], 1);  clear press

%% Setting up components needed for the simulation
% Since we impose no-flow boundary conditions in this example, we restrict
% connections to interior faces only.
N  = double(G.faces.neighbors);
intInx = all(N ~= 0, 2);
N  = N(intInx, :);
% We will be using the two-point flux approximation. First the one-sided
% transmissibilities are computed, then the harmonic average is taken to
% obtain the two-sided transmissibilites.
hT = computeTrans(G, rock);
cf = G.cells.faces(:,1);
nf = G.faces.num;
T  = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]);
T  = T(intInx);
% In setting up the equations, we need discrete forms of the divergence and
% gradient operators, and we represent these as multiplication by sparse
% matrices. In particular, we construct the 'gradient matrix' C as folows:
n = size(N,1);
C = sparse( [(1:n)'; (1:n)'], N, ones(n,1)*[1 -1], n, G.cells.num);
% The discrete gradient and divergence operators are now given by
grad = @(x)-C*x;
div  = @(x)C'*x;

% Additionally, we will need to take the average of cell-based quantities
% in neighboring cells and define the following function
avg = @(x) 0.5 * (x(N(:,1)) + x(N(:,2)));

%% Pressure and well equations:
% The pressure equation (without well contributions) is given by:
%
% $\frac{d}{dt}(\phi\rho)+\nabla\cdot(\rho v)=0, \quad v = -\frac{K}{\mu}\nabla(p-g\rho z)$
%
% In discretized form, this leads to
z = G.cells.centroids(:,3); % z-coordinate of grid cells
pressureEq = @(p, p0, dt) (1/dt) .* (pv(p).*rho(p) - pv(p0).*rho(p0)) ...
    - div( avg(rho(p) ./ mu) .* T .* grad(p - g*rho(p).*z) );

% Wellrates are given as Peaceman well-index times pressure drop
% wc = W(1).cells; % perforation grid cells
% WI = W(1).WI;    % well indices
% dz = W(1).dZ;    % perforation depth relative to well reference depth
% 
% wc2 = W(2).cells; % perforation grid cells
% WI2 = W(2).WI;    % well indices
% dz2 = W(2).dZ;    % perforation depth relative to well reference depth

% wellRates = {...
%    @(p, bhp) W(1).WI .* (rho(p(W(1).cells)) ./ mu) .* (bhp - p( W(1).cells) + g*dz.*rho(p( W(1).cells)));...
%    @(p, bhp) WI2 .* (rho(p(wc2)) ./ mu) .* (bhp - p(wc2) + g*dz2.*rho(p(wc2)))};
wellRates=cell(5,1);
for i=1:5
    wellRates{i}=@(p, bhp) W(i).WI .* (rho(p(W(i).cells)) ./ mu) .* (bhp - p( W(i).cells) + g*W(i).dZ.*rho(p( W(i).cells)));
end
%% Define ADI variables
% The primary variables are grid-cell pressures, well bhp and surface
% rate.

inputv{1}=p_init;
inputv{2}=p_init( W(1).cells);
inputv{3}=p_init( W(2).cells);
inputv{4}=p_init( W(3).cells);
inputv{5}=p_init( W(4).cells);
inputv{6}=p_init( W(5).cells);
inputv{7}=0;
inputv{8}=0;
inputv{9}=0;
inputv{10}=0;
inputv{11}=0;
%ad_values is a vector containing the resulting initialized adi variables
%[output{1}, output{2}, output{3},output{4}, output{5}, output{6},output{7}, output{8},output{9},output{10},output{11}]=
%[     p_ad,    bhp_ad,   bhp_ad2,  bhp_ad3,   bhp_ad4,   bhp_ad5,    qS_ad,    qS_ad2,   qS_ad3,    qS_ad4,    qS_ad5]
 [output{1}, output{2}, output{3},output{4}, output{5}, output{6},output{7}, output{8},output{9},output{10},output{11}]=...
    initVariablesADI(inputv{1},inputv{2},inputv{3},inputv{4},inputv{5},inputv{6},inputv{7},inputv{8},inputv{9},inputv{10},inputv{11});
%return
%[p_ad, bhp_ad, bhp_ad2, qS_ad, qS_ad2] = initVariablesADI(p_init, p_init(wc(1)), p_init(wc2(1)), 0, 0);


%[p_ad,bhp_ad1,bhp_ad2,qS_ad1,qS_ad2]=initVariablesADI(p_init,p_init(wc(1)),p_init(wc(tb)),0,0);
% For convenience, create indices to variables when stacked:
%Ix(1) = pIx, Ix(2)=bhpIx,  Ix(3)=qSIx,Ix(4)=bhpIx2 Ix(5)=qSIx2,
%Ix(6)=bhpIx3, Ix(7)=qSIx3 Ix(8)=bhpIx4, Ix(9)=qSIx4, Ix(10)=bhpIx5,
%Ix(11)= qSIx5 ;
% pIx = 1:G.cells.num; 
% bhpIx  = G.cells.num + 1; qSIx  = G.cells.num + 6;
% bhpIx2 = G.cells.num + 2; qSIx2 = G.cells.num + 7;
% bhpIx3 = G.cells.num + 3; qSIx3 = G.cells.num + 8;
% bhpIx4 = G.cells.num + 4; qSIx4 = G.cells.num + 9;
% bhpIx5 = G.cells.num + 5; qSIx5 = G.cells.num + 10;
Ix{1} = 1:G.cells.num; 
for i=1:10
    Ix{i+1}=G.cells.num+i;
end
    
%% Set up simulation parameters
%
numSteps = 52;
totTime  = 365*day;
dt       = totTime / numSteps;
% Tolerance and maximum number of iterations for the Newton solver.
tol      = 1e-5; 
maxits   = 10;




% Save output in array 'sol'
sol = repmat(struct('time', [], 'pressure', [], 'bhp', [], 'bhp2', [], ...
    'qS', [], 'qS2', []), ...
             [numSteps + 1, 1]);

sol(1).time      = 0;
sol(1).pressure  = double(output{1});
sol(1).bhp       = double(output{2});
sol(1).qS        = double(output{7});
sol(1).bhp2      = double(output{3});
sol(1).qS2       = double(output{8});
sol(1).bhp3      = double(output{4});
sol(1).qS3       = double(output{9});
sol(1).bhp4      = double(output{5});
sol(1).qS4       = double(output{10});
sol(1).bhp5      = double(output{6});
sol(1).qS5       = double(output{11});
% Set up plot
set(0, 'CurrentFigure', f(2));
clf
set(f(2), 'Color', 'w')
subplot(2,1,1),
plotCellData(G, convertTo(p_init, barsa));
title('pressure [bar]', 'EdgeColor', 'w');
colorbar; view(3);
camproj perspective

subplot(2,1,2);
xlim([0, convertTo(totTime,day)]);
title('Surface volume rate [m^3/day]'); 
hold on

temp_loop(dt,totTime,eqs,tol,maxits,output,rhoS,wellRates,Ix,pressureEq,W,f,G,dir)