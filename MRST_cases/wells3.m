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
nx = 10;    ny = 10;    nz = 1;
Dx = 200;   Dy = 200;   Dz = 50;
G = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
G = computeGeometry(G);

%% Setup rock properties
%

% Assume homogeneous/isotropic rock.
permX = 30*milli*darcy;
poro  = 0.3;
rock.perm = repmat(permX, [G.cells.num, 1]);
rock.poro = repmat(poro , [G.cells.num, 1]);

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
p_r  = 00*barsa;

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
p_r   = 200*barsa;
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

W = addWell(W, G, rock, cellInx2, 'Name', 'injector', 'InnerProduct', 'ip_tpf', 'sign', 1);

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

W = addWell(W, G, rock, cellInx, 'Name', 'producer', 'InnerProduct', 'ip_tpf', 'sign', -1);

% Plotting
f = [figure(1), figure(2)];
set(0, 'CurrentFigure', f(1));
clf
plotGrid(G, 'FaceColor', 'g', 'FaceAlpha', .3, 'EdgeColor', 'w');
plotWell(G, W);
axis off;
set(f(1), 'Color', 'w'); 
camproj perspective;
view(3);
pause

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
inputv=cell(5,1);
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
[p_ad, bhp_ad, bhp_ad2, bhp_ad3, bhp_ad4, bhp_ad5,qS_ad, qS_ad2,qS_ad3,qS_ad4,qS_ad5] =...
    initVariablesADI(inputv{1},inputv{2},inputv{3},inputv{4},inputv{5},inputv{6},inputv{7},inputv{8},inputv{9},inputv{10},inputv{11});
%return
%[p_ad, bhp_ad, bhp_ad2, qS_ad, qS_ad2] = initVariablesADI(p_init, p_init(wc(1)), p_init(wc2(1)), 0, 0);


%[p_ad,bhp_ad1,bhp_ad2,qS_ad1,qS_ad2]=initVariablesADI(p_init,p_init(wc(1)),p_init(wc(tb)),0,0);
% For convenience, create indices to variables when stacked:

pIx = 1:G.cells.num; 
bhpIx  = G.cells.num + 1; qSIx  = G.cells.num + 6;
bhpIx2 = G.cells.num + 2; qSIx2 = G.cells.num + 7;
bhpIx3 = G.cells.num + 3; qSIx3 = G.cells.num + 8;
bhpIx4 = G.cells.num + 4; qSIx4 = G.cells.num + 9;
bhpIx5 = G.cells.num + 5; qSIx5 = G.cells.num + 10;
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
sol(1).pressure  = double(p_ad);
sol(1).bhp       = double(bhp_ad);
sol(1).qS        = double(qS_ad);
sol(1).bhp2      = double(bhp_ad2);
sol(1).qS2       = double(qS_ad2);
sol(1).bhp3      = double(bhp_ad3);
sol(1).qS3       = double(qS_ad3);
sol(1).bhp4      = double(bhp_ad4);
sol(1).qS4       = double(qS_ad4);
sol(1).bhp5      = double(bhp_ad5);
sol(1).qS5       = double(qS_ad5);
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


%% Main simulation
% We solve the equations implicitely. At each time step, the equations are assembled and
% the automatic differentiation framework takes care automatically of the computation of
% the Jacobian.

t = 0; 
step = 0;
nDigits = floor(log10(maxits)) + 1;
while t < totTime,
    t = t + dt;
    step = step + 1;
    fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
            step, convertTo(t - dt, day), convertTo(t, day));

    % Newton loop
    resNorm = 1e99;
    p0  = double(p_ad); % Previous step pressure
    nit = 0;
    while (resNorm > tol) && (nit < maxits)
        % Create equations:
        %eqs = cell([3, 1]);

        eqs{1} = pressureEq(p_ad, p0, dt);
        % Add well contributions in perforated cells:
        eqs{1}(W(1).cells) = eqs{1}(W(1).cells) - wellRates{1}(p_ad, bhp_ad);
        eqs{1}(W(2).cells) = eqs{1}(W(2).cells) - wellRates{2}(p_ad, bhp_ad2);
        eqs{1}(W(3).cells) = eqs{1}(W(3).cells) - wellRates{3}(p_ad, bhp_ad3);
        eqs{1}(W(4).cells) = eqs{1}(W(4).cells) - wellRates{4}(p_ad, bhp_ad4);
        eqs{1}(W(5).cells) = eqs{1}(W(5).cells) - wellRates{5}(p_ad, bhp_ad5);

        % Final equation is prescribed bhp at producers (injector)
        eqs{2} = bhp_ad + 100*barsa;
        eqs{3} = bhp_ad2 + 100*barsa;
        eqs{4} = bhp_ad3 + 100*barsa;
        eqs{5} = bhp_ad4 + 100*barsa;
        eqs{6} = bhp_ad5 - 400*barsa;       
        % Sum of wellrates should equal total rate:
        eqs{7} = qS_ad - sum(wellRates{1}(p_ad, bhp_ad))/rhoS;
        eqs{8} = qS_ad2 - sum(wellRates{2}(p_ad, bhp_ad2))/rhoS;
        eqs{9} = qS_ad3 - sum(wellRates{3}(p_ad, bhp_ad3))/rhoS;
        eqs{10} = qS_ad4 - sum(wellRates{4}(p_ad, bhp_ad4))/rhoS;
        eqs{11} = qS_ad5 - sum(wellRates{5}(p_ad, bhp_ad5))/rhoS;
        % Concatenate equations and solve:
        eq  = cat(eqs{:});
        J   = eq.jac{1};  % Jacobian
        res = eq.val;     % residual
        upd = -(J \ res); % Newton update

        % Update variables
        p_ad.val   = p_ad.val   + upd(pIx);
        bhp_ad.val = bhp_ad.val + upd(bhpIx);
        qS_ad.val  = qS_ad.val  + upd(qSIx);
        bhp_ad2.val = bhp_ad2.val + upd(bhpIx2);
        qS_ad2.val  = qS_ad2.val  + upd(qSIx2);
        bhp_ad3.val = bhp_ad3.val + upd(bhpIx3);
        qS_ad3.val  = qS_ad3.val  + upd(qSIx3);
        bhp_ad4.val = bhp_ad4.val + upd(bhpIx4);
        qS_ad4.val  = qS_ad4.val  + upd(qSIx4);
        bhp_ad5.val = bhp_ad5.val + upd(bhpIx5);
        qS_ad5.val  = qS_ad5.val  + upd(qSIx5);
        resNorm = norm(res);
        nit     = nit + 1;
        fprintf('  Iteration %*d:  Res = %.4e\n', nDigits, nit, resNorm);
    end

    if nit > maxits,
        error('Newton solves did not converge')
    else
        sol(step+1).time     = t;
        sol(step+1).pressure = double(p_ad);
        sol(step+1).bhp      = double(bhp_ad);
        sol(step+1).qS       = double(qS_ad);
        sol(step+1).bhp2      = double(bhp_ad2);
        sol(step+1).qS2       = double(qS_ad2);
        sol(step+1).bhp3      = double(bhp_ad3);
        sol(step+1).qS3       = double(qS_ad3);
        sol(step+1).bhp4      = double(bhp_ad4);
        sol(step+1).qS4       = double(qS_ad4);
        sol(step+1).bhp5      = double(bhp_ad5);
        sol(step+1).qS5       = double(qS_ad5);
        % Plot evolution
        set(0, 'CurrentFigure', f(2));
        subplot(2,1,1), cla, %caxis([120, 205])
        plotCellData(G, convertTo(sol(step+1).pressure, barsa), 'EdgeColor', 'w');

        subplot(2,1,2)
        xaxis = convertTo(sol(step+1).time, day);
        y1 = convertTo(-sol(step+1).qS , meter^3/day);
        y2 = convertTo(sol(step+1).qS2 , meter^3/day);
        y3 = convertTo(sol(step+1).qS3 , meter^3/day);
        y4 = convertTo(sol(step+1).qS4 , meter^3/day);
        y5 = convertTo(sol(step+1).qS5 , meter^3/day);
        plot(xaxis, y1,'*', xaxis,y5,'*');
        ylabel('Rate [meter^3/day]','FontSize',14)
    xlabel('Time [days]','FontSize',14)
      legend('W1','W2')
        %plot(); hold on
        drawnow
    end
end
