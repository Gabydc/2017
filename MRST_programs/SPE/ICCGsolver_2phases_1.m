%% Homogeneous quarter five-spot
% In this example, we compare and contrast incompressible single-phase and
% two-phase flow for a homogeneous quarter five-spot. Through the example,
% you will also be introduced to the way well solutions are represented in
% more advanced multiphase simulators based on the AD-OO framework.
close all
clear all
mrstModule add incomp diagnostics
tol_p = 1e-7;
tol_s = 1e-7;
multigrid = false;

linsolve_t = @(J, F) agmg(J, F, 50,tol_s, 2000, 0);


Model_u



hT = computeTrans(G, rock);

% Figure, colormap and contour values
figure('Position',[300 550 1100 650]);
nval = 20;
cval = linspace(0,1,nval+1); cval=.5*cval(1:end-1)+.5*cval(2:end);
colormap(flipud(.5*jet(nval)+.5*ones(nval,3)));
gravity reset off
fluid = initSimpleFluid('mu' , [   1,    1] .* centi*poise     , ...
                        'rho', [1000,  850] .* kilogram/meter^3, ...
                        'n'  , [   2,    2]);

%% Compute solution and plot saturation evolution
% To study the solution, we will plot saturation profiles at four instances
% in time up to slightly beyond water breakthrough; this corresponds
% dimensionless times 0.2 to 0.8 PVI. To get the solution provies as
% accurate as possible, we use the explicit transport solver and M substeps
% to advance the solution 0.2 PVI forward in time. We continue computing
% the solution up to time 1.2 PVI, but do not show snapshots of the
% saturation field for the two last time intervals.

% Compute an initial single-phase pressure solution, from which we estimate
% the final time that corresponds to 1.2 PVI if this flow field remains
% unchanged. With quadratic relperm curves and equal viscosities, the
% multiphase displacement front will propagate at a speed of
% a=1/(2(sqrt(2)-1)) relative to the total velocity.
x  = initState(G,W,500*barsa, [0 1]);
x  = incompTPFA(x, G, hT, fluid, 'wells', W);
xr  = incompTPFA(x, G, hT, fluid, 'wells', W);
T  = 1.2*sum(pv)/x.wellSol(1).flux;
a  = 1/(2*(sqrt(2)-1));

% Compute time-of-flight for the single-phase flow field and record the
% corresponding breakthrough time in the producer.
tau = computeTimeOfFlight(x, G, rock,  'wells', W);
tbf = tau(W(2).cells,1);

% Initialize number of time intervals, cell array to hold well solutions,
% and array to hold the oil in place
[N,M]    = deal(10,10);
step = T / N*M;

wellSols = cell(N*M+1,1);  wellSols{1} = getWellSol(W, x, fluid);
oip      = zeros(N*M+1,1); oip(1) = sum(x.s(:,2).*pv);

wellSolsr = cell(N*M+1,1);  wellSolsr{1} = getWellSol(W, xr, fluid);
oipr      = zeros(N*M+1,1); oipr(1) = sum(xr.s(:,2).*pv);

Changing_w

dir = '/mnt/sda2/cortes/Programs/2017/MRST_programs/SPE/Results1/';
for n=1:N
    
    fprintf(1,'Main step %d: ',n);
    for m=1:M
step = (n-1)*N+m
        W=W1{step};
p0 = x.pressure;
np = size(p0);
for i=1:numel(W)
    p0(np+i) = 0;
end
if (multigrid == true)
linsolve_p = @(S, h) agmg(S, h,  1, tol_p, 1000, 0);
%agmg(A,b,icg,tol,maxit,verbose,x0,ijob)
else
solver = ICCGSolverAD('tolerance', tol_p,'maxIterations', 1000,'x0',p0,...
    'dir',[],'W', W);
linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
end
        
        if (multigrid == true)
            psolve = @(x,W)incompTPFA_amg(x, G, hT, fluid, 'wells', W, 'LinSolve', linsolve_p);
        else
            psolve = @(x,W)incompTPFA_Def(x, G, hT, fluid, 'wells', W,'LinSolve', linsolve_p);
        end
        tic
        [x] = psolve(x,W);
        tsol=toc;
        tic
       xr  = incompTPFA(xr, G, hT, fluid, 'wells', W);
       sol=toc;
         sts = 1;
       for si= 1:sts
       x  = implicitTransport(x, G, T/(N*M*si), rock, fluid, 'wells', W, 'LinSolve', linsolve_t);
       xr  = implicitTransport(xr, G, T/(N*M*si), rock, fluid, 'wells', W, 'LinSolve', linsolve_t);
       end
      
       wellSols{(n-1)*M+m+1} = getWellSol(W, x, fluid);
       oip((n-1)*M+m+1) = sum(x.s(:,2).*pv);
       if x.s(W(2).cells,1)<eps
           W(2).bt = (n-1)*M+m+1;
       end
       fprintf(1,'%d, ',m);
         fprintf(1,' ');
         
         
                  wellSolsr{(n-1)*M+m+1} = getWellSol(W, xr, fluid);
       oipr((n-1)*M+m+1) = sum(xr.s(:,2).*pv);
       if xr.s(W(2).cells,1)<eps
           Wr(2).bt = (n-1)*M+m+1;
       end
       fprintf(1,'%d, ',m);
    
    fprintf(1,'\n');
       
    end
  
   
    
end
%% Plotting with GUI from ad-core
mrstModule add ad-core
plotWellSols(wellSols,cumsum(dt))







