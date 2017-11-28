%% Simulate the SPE10 base case
% The simulation uses a mimetic pressure solver and an implicit transport
% solver.

clear, close all hidden

%%
use_mimetic = false;
use_reorder = false;
use_ICCG = false;
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

cartDims = [60, 220, numel(layers)];
rock     = SPE10_rock(layers);

rock.perm = convertFrom(rock.perm, milli*darcy);

is_pos             = rock.poro > 0;
rock.poro(~is_pos) = min(rock.poro(is_pos));

physDims = cartDims .* [20, 10, 2]*ft;

G = computeGeometry(cartGrid(cartDims, physDims));

if use_mimetic,
   try
      require mimetic
   catch
      mrstModule add mimetic
   end

   S = computeMimeticIP(G, rock);
else
   T = computeTrans(G, rock);
end

%%
fluid = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                        'rho', [1014, 859]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);

%%
% Set Comp_i=[0,0] in producers to counter X-flow effects...
%
if use_mimetic,
   well_ip = 'ip_simple';
else
   well_ip = 'ip_tpf';
end
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

%% Change wells 
% Changing_w

%%
x         = initResSol (G, 0);
x.wellSol = initWellSol(W, 0);

%%
linsolve_p = @(S, h) agmg(S, h,  1, 5.0e-11, 1000, 0);
linsolve_t = @(J, F) agmg(J, F, 50, 5.0e-11, 2000, 0);

if use_mimetic,
   psolve = @(x) ...
      incompMimetic(x, G, S, fluid, 'wells', W, 'LinSolve', linsolve_p);
else
   psolve = @(x) ...
      incompTPFA(x, G, T, fluid, 'wells', W, 'LinSolve', linsolve_p);
end

if ~use_reorder,
   tsolve = @(x, dt) ...
      implicitTransport(x, G, dt, rock, fluid, 'wells', W, ...
                        'LinSolve', linsolve_t,'verbose',true);
else
   mu          = fluid.properties();
   fluid.param = struct('viscw', mu(1), 'visco', mu(2),       ...
                        'srw', 0, 'sro', 0, 'nw', 1, 'no', 1, ...
                        'satnum', ones([G.cells.num, 1], 'int32'));

   tsolve = @(x, dt) ...
      implicitTransportReorder(x, G, dt, rock, fluid, 'wells', W);
end

%%
DT    = 100*day;
nstep =  10;      % 2000 days

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

for k = 1 : nstep,
   t0 = tic; x = psolve(x);     dt = toc(t0);
   fprintf('[%02d]: Pressure:  %12.5f [s]\n', k, dt);

   t0 = tic; x = tsolve(x, DT); dt = toc(t0);
   fprintf('[%02d]: Transport: %12.5f [s]\n', k, dt);

   t = k * DT;

   [wres{:}] = prodCurves(W, x, fluid);
   Prod      = append_wres(Prod, t, wres{:});
end

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
               % subplotcbspe(nf1,clim,k,Np,G,nz,time,Press)





