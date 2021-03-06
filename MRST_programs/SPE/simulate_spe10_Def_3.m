%% Simulate the SPE10 base case
% The simulation uses a mimetic pressure solver and an implicit transport
% solver.

clear, close all hidden

%%
dir='/mnt/sda2/cortes/Results/2017/Report/SPE10/training/11_26/ex1/';
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
layers = 1 : 85;
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
                 'Val', 7.6590e7, 'Radius', 0.125*meter, ...
                 'Name', 'I1', 'Comp_i', [1, 0]);



%%
x         = initResSol (G, 20000*psia);
x.wellSol = initWellSol(W, 0);

%%
linsolve_p = @(S, h) agmg(S, h,  1, 5.0e-11, 1000, 0);
psolve = @(x) ...
    incompTPFA(x, G, T, fluid, 'wells', W, 'LinSolve', linsolve_p);

linsolve_t = @(J, F) agmg(J, F, 50, 5.0e-7, 2000, 0);
tsolve = @(x, dt) ...
    implicitTransport(x, G, dt, rock, fluid, 'wells', W, ...
    'LinSolve', linsolve_t,'verbose',true);

%%
DT    = 100*day;
nstep =  40;
%% Change wells

% Pressure in injector and producers
P = 275;
I = 700;

folder=[ 'SPE10_' num2str(numel(layers))  'DT_' num2str(DT/day) 'step_' num2str(nstep) 'P_1'];
mkdir([dir], folder)
dir1 = [dir folder '/'];
if use_ICCG
    def = 0;
    training = 0;
    folder=['ICCG' ];
else
    def = 1;
    pod = 1;
    dv = 10;
    training = 1;
    folder=['DICCG' ];
end

mkdir([dir1], folder)
dir2 = [dir1 folder '/'];




% Number of time steps with same pressure
tch =2;
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
tol = 5.0e-8;
maxIter = 1500;
%%
nf=0;
nf1=nf+1;
n=nx*ny*nz;

if(use_DICCG)
    files=['Pressure'];
    filename=[dir1 files ];
    load(filename)
    
    np = size(Pressure,2);
    dpod = [np-dv+1:np];
    [U,S]=PODbasis(Pressure);
    Z=U(:,dpod);
    nf = nf+1;
    file{nf} = ['eig_pod'];
    f(nf) = figure(nf);
    plot(log((diag(S))),'*r');
    ylabel('log(Value) ','FontSize',16)
    xlabel('Eigenvalue','FontSize',16)
    axis('tight');
    nf = nf + 1;
    
    file{nf} = ['eig_vect' ];
    f(nf) = figure(nf);
    plot(Z);
    ylabel('vectors ','FontSize',16)
    xlabel('eigenvectors','FontSize',16)
    axis('tight')
    
    for i=1:numel(W)
        Z(n+i,1)=0;
    end
end


[wres{:}]        = prodCurves(W, x, fluid);
Pressure         = zeros(n,nstep);
fluxes     = zeros(length(x.flux),nstep);
Sat       = zeros(n,nstep);
facePressure= zeros(length(x.flux),nstep);

load I_p
dI_p = (I_p(20) - I_p(1))/nstep;
for k = 1 : nstep,
    
    PI_1 = (I_p(1)+dI_p*k)*0.5;
    W(5).val = PI_1;
    
    % Set initial condition
    p0 = x.pressure;
    
    for i=1:numel(W)
        p0(n+i) = 0;
    end
    
    if use_DICCG,
        diccg=1
        W = W0{k};
         W(5).val = PI_1;
        solver = DICCGSolverAD('tolerance', tol,'maxIterations',  maxIter,'Z',Z,'x0',p0,'W', W);
        linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
        psolve = @(x) incompTPFA_Def(x, G, T, fluid, 'wells', W,'LinSolve', linsolve_p);
        t0 = tic;
        [x,preport(k)]= psolve(x);
        dt = toc(t0);
    else if use_ICCG
            iccg =1
            W = W1{k};
            W(5).val = PI_1;
            solver = ICCGSolverAD('tolerance', tol,'maxIterations',  maxIter,'x0',p0,'W', W);
            linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
            psolve = @(x) incompTPFA_Def(x, G, T, fluid, 'wells', W,'LinSolve', linsolve_p);
            t0 = tic;
            [x,preport(k)]= psolve(x);
            dt = toc(t0);
        else
            mg = 1
            linsolve_p = @(S, h) agmg(S, h,  1, tol, maxIter, 0);
            psolve = @(x) incompTPFA(x, G, T, fluid, 'wells', W, 'LinSolve', linsolve_p);
            t0 = tic;
            x = psolve(x);
            dt = toc(t0);
        end
    end
    fprintf('[%02d]: Pressure:  %12.5f [s]\n', k, dt);
    
    t0 = tic;
    % for i= 1:10
    x = tsolve(x, DT);
    % end
    dt = toc(t0);
    fprintf('[%02d]: Transport: %12.5f [s]\n', k, dt);
    
    t = k * DT;
    
    [wres{:}]        = prodCurves(W, x, fluid);
    Pressure(:,k)    = x.pressure;
    fluxes(:,k)      = x.flux;
    Sat(:,k)         = x.s;
    [wellSol{k}]     = x.wellSol;
    facePressure(:,k)= x.facePressure;
    Prod             = append_wres(Prod, t, wres{:});
    Prod1{k}        = Prod;
    for i=1:5
        pw1(i,k)=x(1).wellSol(i).pressure;
    end
end
%%
ex=1;
pod =1;
dv = 20;
np = size(Pressure,2);
dpod = [np-dv+1:np];
filetx = ['results.txt'];
saveits(dir1,filetx,def,pod,dpod,k,dv,preport)


clear figure
filews=['workspace'];
filename=[dir2 filews];
save(filename)



%%
nf = nf + 1;
figure(nf);
file{nf} = ['Pressure'];
for i=1:4
    plot(pw1(i,:)/barsa,'*-','color', [0.1*i 0.5 0.6])
    hold on
end

%figure; plot(wi/barsa)

%%
nf = nf + 1;
figure(nf);
file{nf} = ['Sat'];
plotCellData(G, x.s(:,1), 'EdgeColor', 'k', ...
    'EdgeAlpha', 0.050, 'FaceAlpha', 0.375)
view(3), colorbar, axis tight off

nf = nf + 1;
figure(nf);
file{nf} = ['Sat_1'];
plotCellData(G, x.s(:,1), find(x.s > 0.5), 'EdgeColor', 'k', ...
    'EdgeAlpha', 0.050, 'FaceAlpha', 0.375)
view(3), colorbar, axis tight off

nf = nf + 1;
figure(nf);
file{nf} = ['Tot_prod_rate'];
plot(convertTo(Prod.t, day), convertTo(Prod.vpt(:,1:end-1), meter^3/day))
legend({ W(1:end-1).name }, 'Location', 'Best')
xlabel('Time [d]'), ylabel('Total Production Rate [m^3/d]')

nf = nf + 1;
figure(nf);
file{nf} = ['Oil_prod_rate'];
plot(convertTo(Prod.t, day), convertTo(Prod.opr(:,1:end-1), meter^3/day))
legend({ W(1:end-1).name }, 'Location', 'Best')
xlabel('Time [d]'), ylabel('Oil Production Rate [m^3/d]')

nf = nf + 1;
figure(nf);
file{nf} = ['Water_prod_rate'];
plot(convertTo(Prod.t, day), convertTo(Prod.wpr(:,1:end-1), meter^3/day))
legend({ W(1:end-1).name }, 'Location', 'Best')
xlabel('Time [d]'), ylabel('Water Production Rate [m^3/d]')

nf = nf + 1;
figure(nf);
file{nf} = ['Well_wat_cut'];
plot(convertTo(Prod.t,day), Prod.wc(:,1:end-1))
legend({ W(1:end-1).name }, 'Location', 'Best')
xlabel('Time [d]'), ylabel('Well Water Cut')

%%
Press = x.pressure;
Sat_w = x.s;

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
nf = nf + 1;
figure(nf);
%title('Pressure [bars]');
file{nf} = ['Pressure'];
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




%%


%%

for i = nf1 : nf
    f(i) = figure(i);
    savefigures(f(i), file{i}, dir2)
end
clear f
%%


% files=['Pressure'];
% filename=[dir1 files ];
% save(filename, 'Pressure')



