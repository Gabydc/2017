%% Simulate the SPE10 base case
% The simulation uses a mimetic pressure solver and an implicit transport
% solver.

clear, close all hidden

%%

use_ICCG = true;
use_DICCG = false;

spe10_data  = fullfile(fileparts(mfilename('fullpath')), ...
                       '..', 'spe10_rock.mat');
if ~exist(spe10_data, 'file'),
   if ~make_spe10_data,
      error(['Failed to establish on-disk representation of ', ...
             'SPE10 rock data']);
   end
end

%%
layers = 1 : 35;
[nx, ny, nz] = deal(60, 220, numel(layers));
cartDims = [nx, ny,nz];
rock     = SPE10_rock(layers);

rock.perm = convertFrom(rock.perm, milli*darcy);

is_pos             = rock.poro > 0;
rock.poro(~is_pos) = min(rock.poro(is_pos));

physDims = cartDims .* [20, 10, 2]*ft;

G = computeGeometry(cartGrid(cartDims, physDims));

 T = computeTrans(G, rock);


%%
fluid = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                        'rho', [1014, 859]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);

%%
% Set Comp_i=[0,0] in producers to counter X-flow effects...
%
well_ip = 'ip_tpf';
W = verticalWell([], G, rock,  1,   1, [], 'Type', 'bhp', ...
                 'InnerProduct', well_ip, ...
                 'Val', 4000*psia, 'Radius', 0.125*meter, ...
                 'Name', 'P1', 'Comp_i', [0, 0]);

W = verticalWell(W , G, rock, 60,   1, [], 'Type', 'bhp', ...
                 'InnerProduct', well_ip, ...
                 'Val', 4000*psia, 'Radius', 0.125*meter, ...
                 'Name', 'P2', 'Comp_i', [0, 0]);

W = verticalWell(W , G, rock, 60, 220, [], 'Type', 'bhp', ...
                 'InnerProduct', well_ip, ...
                 'Val', 4000*psia, 'Radius', 0.125*meter, ...
                 'Name', 'P3', 'Comp_i', [0, 0]);

W = verticalWell(W , G, rock,  1, 220, [], 'Type', 'bhp', ...
                 'InnerProduct', well_ip, ...
                 'Val', 4000*psia, 'Radius', 0.125*meter, ...
                 'Name', 'P4', 'Comp_i', [0, 0]);

W = verticalWell(W , G, rock, 30, 110, [], 'Type', 'bhp',   ...
                 'InnerProduct', well_ip, ...
                 'Val', 40000*psia, 'Radius', 0.125*meter, ...
                 'Name', 'I1', 'Comp_i', [1, 0]);



%%
x         = initResSol (G, 0);
x.wellSol = initWellSol(W, 0);

%%
linsolve_p = @(S, h) agmg(S, h,  1, 5.0e-11, 1000, 0);
psolve = @(x) ...
      incompTPFA(x, G, T, fluid, 'wells', W, 'LinSolve', linsolve_p);
  
linsolve_t = @(J, F) agmg(J, F, 50, 5.0e-11, 2000, 0);
   tsolve = @(x, dt) ...
      implicitTransport(x, G, dt, rock, fluid, 'wells', W, ...
                        'LinSolve', linsolve_t,'verbose',true);

%%
DT    = 10*day;
nstep =  60;      
%% Change wells 

% Pressure in injector and producers
P = 4000;
I = 40000;
% Number of time steps with same pressure
 tch =6;

 Changing_w
 %for i=1:6; wi(i)=W1{i}(1).val; end 
 %figure; plot(wi/barsa)
%break
Prod = struct('t'  , []                  , ...
              'vpt', zeros([0, numel(W)]), ...
              'opr', zeros([0, numel(W)]), ...
              'wpr', zeros([0, numel(W)]), ...
              'wc' , zeros([0, numel(W)]));

append_wres = @(x, t, vpt, opr, wpr, wc) ...
   struct('t'  , [x.t  ; t                  ], ...
          'vpt', [x.vpt; reshape(vpt, 1, [])], ...
          'opr', [x.opr; reshape(opr, 1, [])], ...
          'wpr', [x.wpr; reshape(wpr, 1, [])], ...
          'wc' , [x.wc ; reshape(wc , 1, [])]);

wres = cell([1, 4]);
pw1=zeros(5,nstep);
tol = 5.0e-11;
maxIter = 1000;

for k = 1 : nstep,
    W = W1{k};
    % Set initial condition
    p0 = x.pressure;
n=nx*ny*nz;
for i=1:numel(W)
    p0(n+i) = 0;
end
    if use_ICCG
        iccg =1
        solver = ICCGSolverAD('tolerance', tol,'maxIterations',  maxIter,'x0',p0,'W', W);
        linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
        psolve = @(x) incompTPFA_Def(x, G, T, fluid, 'wells', W,'LinSolve', linsolve_p);
        t0 = tic;
        [x,preport(k)]= psolve(x);
        dt = toc(t0);
    else
        mg = 1
        linsolve_p = @(S, h) agmg(S, h,  1, tol, 1000, 0);
        psolve = @(x) incompTPFA(x, G, T, fluid, 'wells', W, 'LinSolve', linsolve_p);
        t0 = tic;
        x = psolve(x);
        dt = toc(t0);
    end
   fprintf('[%02d]: Pressure:  %12.5f [s]\n', k, dt);

   t0 = tic; 
   x = tsolve(x, DT); 
   dt = toc(t0);
   fprintf('[%02d]: Transport: %12.5f [s]\n', k, dt);

   t = k * DT;

   [wres{:}] = prodCurves(W, x, fluid);
   [xsol{k}] = x;
   Prod      = append_wres(Prod, t, wres{:});
for i=1:5
     pw1(i,k)=x(1).wellSol(i).pressure;
end
end
%%
figure
for i=1:4
    plot(pw1(i,:)/barsa,'*-','color', [0.1*i 0.5 0.6])  
    hold on
end
 %figure; plot(wi/barsa)

%%
figure
plotCellData(G, x.s(:,1), 'EdgeColor', 'k', ...
             'EdgeAlpha', 0.050, 'FaceAlpha', 0.375)
view(3), colorbar, axis tight off

figure
plotCellData(G, x.s(:,1), find(x.s > 0.5), 'EdgeColor', 'k', ...
             'EdgeAlpha', 0.050, 'FaceAlpha', 0.375)
view(3), colorbar, axis tight off

figure
plot(convertTo(Prod.t, day), convertTo(Prod.vpt(:,1:end-1), meter^3/day))
legend({ W(1:end-1).name }, 'Location', 'Best')
xlabel('Time [d]'), ylabel('Total Production Rate [m^3/d]')

figure
plot(convertTo(Prod.t, day), convertTo(Prod.opr(:,1:end-1), meter^3/day))
legend({ W(1:end-1).name }, 'Location', 'Best')
xlabel('Time [d]'), ylabel('Oil Production Rate [m^3/d]')

figure
plot(convertTo(Prod.t, day), convertTo(Prod.wpr(:,1:end-1), meter^3/day))
legend({ W(1:end-1).name }, 'Location', 'Best')
xlabel('Time [d]'), ylabel('Water Production Rate [m^3/d]')

figure
plot(convertTo(Prod.t,day), Prod.wc(:,1:end-1))
legend({ W(1:end-1).name }, 'Location', 'Best')
xlabel('Time [d]'), ylabel('Well Water Cut')

%%
Press = x.pressure;
 Sat_w = x.s;
 nf=1000;
               pmin=0.9*min(Press);
                pmax=1.2*max(Press);
                
                plotNo = 0;
                Np = 3;
                time = (0:DT:DT*k)/day;
                px = [0.01 0.26 0.51 0.76];
                
                nf = nf + 1;
                figure(nf);
                %title('Water Saturation');
                file{nf} = ['Water_saturation'];
                clim = [0 1];
                nz = numel(layers);
                plotCellData(G, Sat_w,'LineStyle','none', 'EdgeColor', 'k', ...
             'EdgeAlpha', 0.050, 'FaceAlpha', 0.375);
                for i = 1 : 4
                    plotWell(G, W(i),'color', 'r');
                end
                plotWell(G, W(5), 'color', 'b');
                view(3)
                axis equal off 
                colorbar('south')
               % subplotcbspe(nf,clim,k,Np,G,nz,time,Sat_w)
                % subplotcbwbt(nf,clim,wbt,Np,G,nz,time,Sat_w)
                nf1 = nf + 1;
                figure(nf1);
                %title('Pressure [bars]');
                file{nf1} = ['Pressure'];
                clim = [pmin pmax];
                plotCellData(G, Press,'LineStyle','none');
                for i = 1 : 4
                    plotWell(G, W(i),'color', 'r');
                end
                plotWell(G, W(5), 'color', 'b');
                view(3)
                axis equal off 
                colorbar('south')
               % subplotcbspe(nf1,clim,k,Np,G,nz,time,Press)

clear figure
save workspace


