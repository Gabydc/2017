%% Simulate the SPE10 base case
% The simulation uses a mimetic pressure solver and an implicit transport
% solver.

clear, close all hidden

%%
use_mimetic = false;
use_reorder = false;
spe10_data  = fullfile(fileparts(mfilename('fullpath')), ...
                       '..', 'spe10_rock.mat');
if ~exist(spe10_data, 'file'),
   if ~make_spe10_data,
      error(['Failed to establish on-disk representation of ', ...
             'SPE10 rock data']);
   end
end

%%
layers = 1 : 85;

cartDims = [60, 220, numel(layers)];
rock     = SPE10_rock(layers);

rock.perm = convertFrom(rock.perm, milli*darcy);

is_pos             = rock.poro > 0;
rock.poro(~is_pos) = min(rock.poro(is_pos));

physDims = cartDims .* [6, 3, 0.6]*meter;

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

%load I_p
if use_mimetic,
   well_ip = 'ip_simple';
else
   well_ip = 'ip_tpf';
end
W = verticalWell([], G, rock,  1,   1, [], 'Type', 'bhp', ...
                 'InnerProduct', well_ip, ...
                 'Val', 275*barsa, 'Radius', 0.125*meter, ...
                 'Name', 'P1', 'Comp_i', [0, 0]);

W = verticalWell(W , G, rock, 60,   1, [], 'Type', 'bhp', ...
                 'InnerProduct', well_ip, ...
                 'Val', 275*barsa, 'Radius', 0.125*meter, ...
                 'Name', 'P2', 'Comp_i', [0, 0]);

W = verticalWell(W , G, rock, 60, 220, [], 'Type', 'bhp', ...
                 'InnerProduct', well_ip, ...
                 'Val', 275*barsa, 'Radius', 0.125*meter, ...
                 'Name', 'P3', 'Comp_i', [0, 0]);

W = verticalWell(W , G, rock,  1, 220, [], 'Type', 'bhp', ...
                 'InnerProduct', well_ip, ...
                 'Val', 275*barsa, 'Radius', 0.125*meter, ...
                 'Name', 'P4', 'Comp_i', [0, 0]);

W = verticalWell(W , G, rock, 30, 110, [], 'Type', 'bhp',   ...
                 'InnerProduct', well_ip, ...
                 'Val',  7.6590e7, 'Radius', 0.125*meter, ...
                 'Name', 'I1', 'Comp_i', [1, 0]);

%%
x         = initResSol (G, 0);
x.wellSol = initWellSol(W, 0);

%%
linsolve_p = @(S, h) agmg(S, h,  1, 5.0e-11, 1000, 0);
linsolve_t = @(J, F) agmg(J, F, 50, 5.0e-11, 1000, 0);

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
                        'LinSolve', linsolve_t);
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
nstep =  20;      % 2000 days

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

load I_p
dI_p = (I_p(20) - I_p(1))/nstep;
for k = 1 : nstep,
    
    PI_1 = (I_p(1)+dI_p*k)*0.5;
    W(5).val = PI_1;
    
    %pause
    linsolve_p = @(S, h) agmg(S, h,  1, 5.0e-11, 1000, 0);
    psolve = @(x) incompTPFA(x, G, T, fluid, 'wells', W, 'LinSolve', linsolve_p);
    t0 = tic;
    x = psolve(x);
    dt = toc(t0);
    
  % t0 = tic; x = psolve(x);     dt = toc(t0);
   fprintf('[%02d]: Pressure:  %12.5f [s]\n', k, dt);

   t0 = tic; x = tsolve(x, DT); dt = toc(t0);
   fprintf('[%02d]: Transport: %12.5f [s]\n', k, dt);

   t = k * DT;

   [wres{:}] = prodCurves(W, x, fluid);
   Prod      = append_wres(Prod, t, wres{:});
   [I_p2(k)]     = x.wellSol(5).pressure;
  % PI_1 =I_p1(k)*0.8
end
%%
figure
plotCellData(G, x.s(:,1),'LineStyle','none')
view(3), colorbar, axis equal tight off
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
