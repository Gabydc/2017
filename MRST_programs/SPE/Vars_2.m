%% Open SPE 10 model
spe10_data  = fullfile(fileparts(mfilename('fullpath')), ...
    '..', 'spe10_rock.mat');
if ~exist(spe10_data, 'file'),
    if ~make_spe10_data,
        error(['Failed to establish on-disk representation of ', ...
            'SPE10 rock data']);
    end
end

%% Initial values
units_w = barsa;
%Pressure of the wells
P = 275;
I = 500;
 
%Pressure of the reservoir
p_0 = 0;

%% Time
% Time steps
DT    = 50*day;
nstep =  10;

%% Model
layers = 1 : 35;
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
    'rho', [1014, 859]*kilogram/meter^3, ...
    'n'  , [   2,   2]);


%% Define Solver
use_ICCG   = true;
use_DICCG  = false;
use_POD    = true;
plot_sol   = true; 
save_res   = true;
training   = true;
% Solvers variables
tol = 5.0e-8;
maxIter = 1500;

% Deflation parameters
dv = 15;
dpod = [nstep-dv+1:nstep];
last       = false;
%% Create directories to save results
dir='/mnt/sda2/cortes/Results/2017/Report/SPE10/training/12_03/ex1/';
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

W = verticalWell(W , G, rock, ceil(nx/2), ceil(ny/2), [], 'Type', 'bhp',   ...
                 'InnerProduct', well_ip, ...
                 'Val', I*barsa, 'Radius', 0.125*meter, ...
                 'Name', 'I1', 'Comp_i', [1, 0]);
%% Changing wells parameters

% Number of time steps with same pressure
tch =2;
Changing_w
