%% Simulate the SPE10 base case
% The simulation uses a mimetic pressure solver and an implicit transport
% solver.

clear, close all hidden



%% Initial values
units_P = barsa;
%Pressure of the wells
P = 270;
I = 5000;
%I = 50;
units_I = stb/day;
%units_I = barsa;
%Pressure of the reservoir
%p_0 = 300;

%% Time
% Time steps
DT    = 50*day;
nstep =  5;
%% Model
layers = 1 : 85;
[nx, ny, nz] = deal(60, 220, numel(layers));
cartDims = [nx, ny,nz];
physDims = cartDims .* [20, 10, 2]*ft;
N = nx*ny*nz;

% Construct the model
rock     = SPE10_rock(layers);
rock.perm = convertFrom(rock.perm, milli*darcy);
is_pos             = rock.poro > 0;
rock.poro(~is_pos) = min(rock.poro(is_pos));
G = computeGeometry(cartGrid(cartDims, physDims));
T = computeTrans(G, rock);

%% Fluids properties
fluid = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
    'rho', [1000, 700]*kilogram/meter^3, ...
    'n'  , [   2,   2]);

%% Compute wells 
% Set Comp_i=[0,0] in producers to counter X-flow effects...

well_ip = 'ip_tpf';
W = verticalWell([], G, rock,  1,   1, [], 'Type', 'bhp', ...
                 'InnerProduct', well_ip, ...
                 'Val', P*barsa, 'Radius', 0.125*meter, ...
                 'Name', 'P1', 'Comp_i', [0, 0]);

W = verticalWell(W , G, rock, nx,   1, [], 'Type', 'bhp', ...
                 'InnerProduct', well_ip, ...
                 'Val', P*barsa, 'Radius', 0.125*meter, ...
                 'Name', 'P2', 'Comp_i', [0, 0]);

W = verticalWell(W , G, rock, nx, ny, [], 'Type', 'bhp', ...
                 'InnerProduct', well_ip, ...
                 'Val', P*barsa, 'Radius', 0.125*meter, ...
                 'Name', 'P3', 'Comp_i', [0, 0]);

W = verticalWell(W , G, rock,  1, ny, [], 'Type', 'bhp', ...
                 'InnerProduct', well_ip, ...
                 'Val', P*barsa, 'Radius', 0.125*meter, ...
                 'Name', 'P4', 'Comp_i', [0, 0]);

W = verticalWell(W , G, rock, ceil(nx/2), ceil(ny/2), [], 'Type', 'rate',   ...
                 'InnerProduct', well_ip, ...
                 'Val', I*units_I, 'Radius', 0.125*meter, ...
                 'Name', 'I1', 'Comp_i', [1, 0]);
             %% Changing wells parameters

% Number of time steps with same pressure
tch =2;
Changing_w
%%
x         = initResSol (G, 0);
x.wellSol = initWellSol(W, 0);

%%
linsolve_p = @(S, h) agmg(S, h,  1, 5.0e-11, 1000, 0);
linsolve_t = @(J, F) agmg(J, F, 50, 5.0e-11, 2000, 0);


   psolve = @(x) ...
      incompTPFA(x, G, T, fluid, 'wells', W, 'LinSolve', linsolve_p);


   tsolve = @(x, dt) ...
      implicitTransport(x, G, dt, rock, fluid, 'wells', W, ...
                        'LinSolve', linsolve_t, 'verbose',true);


%% Define Solver
use_ICCG   = true;
use_DICCG  = false;
use_POD    = true;
plot_sol   = false; 
save_res   = true;
training   = true;
% Solvers variables
tol = 5.0e-8;
maxIter = 2000;

% Deflation parameters
dv = 20;
dpod = [nstep-dv+1:nstep];
last       = false;

%% Create directories to save results
dir='../../12_06/';
folder=[ 'SPE10_' num2str(numel(layers))  'DT_' num2str(DT/day) 'step_' num2str(nstep) 'P_1'];
mkdir([dir], folder)
dir1 = [dir folder '/'];
if use_ICCG
    folder=['ICCG' ];
else
    folder=['DICCG' ];
end
mkdir([dir1], folder)
dir2 = [dir1 folder '/'];
