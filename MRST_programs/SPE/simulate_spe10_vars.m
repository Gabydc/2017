%% Simulate the SPE10 base case
% The simulation uses a mimetic pressure solver and an implicit transport
% solver.

clear, close all hidden
%% Define all the variables

Vars_1


%%
x         = initResSol (G, p_0*barsa);
x.wellSol = initWellSol(W, 0);

%%
linsolve_p = @(S, h) agmg(S, h,  1, tol, maxIter, 0);
psolve = @(x) ...
    incompTPFA(x, G, T, fluid, 'wells', W, 'LinSolve', linsolve_p);

linsolve_t = @(J, F) agmg(J, F, 50, tol, maxIter, 0);
tsolve = @(x, dt) ...
    implicitTransport(x, G, dt, rock, fluid, 'wells', W, ...
    'LinSolve', linsolve_t,'verbose',true);

%% Change wells


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

%%

if(use_DICCG)
    use_DICCG
    files=['Pressure'];
    filename=[dir1 files ];
    load(filename)
    np = size(Pressure,2);
    dpod = [np-dv+1:np];
    [U,S]=PODbasis(Pressure);
    Z=U(:,dpod);
    
    for i=1:numel(W)
        Z(N+i,1)=0;
    end
end


[wres{:}]        = prodCurves(W, x, fluid);
Pressure         = zeros(N,nstep);
fluxes           = zeros(length(x.flux),nstep);
Sat              = zeros(N,nstep);
facePressure     = zeros(length(x.flux),nstep);

%load I_p
%dI_p = (I_p(20) - I_p(1))/nstep;
for k = 1 : nstep,
    
    %PI_1 = I_p(1)+dI_p*k;
    %W(5).val = PI_1;
    
    % Set initial condition
    p0 = x.pressure;
    
    for i=1:numel(W)
        p0(N+i) = 0;
    end
    
    if(use_DICCG)
        use_DICCG
        W = W0{k};
        %  W(5).val = PI_1;
        solver = DICCGSolverAD('tolerance', tol,'maxIterations',  maxIter,'Z',Z,'x0',p0,'W', W);
        linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
        psolve = @(x) incompTPFA_Def(x, G, T, fluid, 'wells', W,'LinSolve', linsolve_p);
        t0 = tic;
        [x,preport(k)]= psolve(x);
        dt = toc(t0);
        
    else if(use_ICCG)
            use_ICCG
            W = W1{k};
            %   W(5).val = PI_1;
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
    x = tsolve(x, DT);
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
    Prod1{k}         = Prod;
    for i=1:5
        pw1(i,k)=x(1).wellSol(i).pressure;
    end
end


%% Save results
if save_res
    if(~training)
        for i=1:ts
            its(i,1)=preport(1,i).iter;
        end
        ttits = sum(its);
        save([dir1  'ttits_T.mat'],'ttits_T')
    else
        filetx = ['resultst.txt'];
        saveits(dir1,filetx,use_ICCG,use_DICCG,use_POD,dpod,k,dv,preport,last)
    end
    filews=['workspace'];
    filename=[dir2 filews];
    save(filename)
    %  clearvars -except Pressure dir1 plot_sol
    if training
        filews=['Pressure'];
        filename=[dir1 filews];
        save(filename,'Pressure')
    end
end
%%
if plot_sol
    Plot_1
end



