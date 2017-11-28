%% Homogeneous quarter five-spot
% In this example, we compare and contrast incompressible single-phase and
% two-phase flow for a homogeneous quarter five-spot. Through the example,
% you will also be introduced to the way well solutions are represented in
% more advanced multiphase simulators based on the AD-OO framework.
close all
clear all
mrstModule add incomp diagnostics
tol_p = 1e-8;
multigrid = false;
linsolve_t = @(J, F) agmg(J, F, 50,tol_s, 2000, 0);




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
x  = initState(G,W,100*barsa, [0 1]);
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
[N,M]    = deal(6,10);
wellSols = cell(N*M+1,1);  wellSols{1} = getWellSol(W, x, fluid);
oip      = zeros(N*M+1,1); oip(1) = sum(x.s(:,2).*pv);

wellSolsr = cell(N*M+1,1);  wellSolsr{1} = getWellSol(W, xr, fluid);
oipr      = zeros(N*M+1,1); oipr(1) = sum(xr.s(:,2).*pv);


dir = '/mnt/sda2/cortes/Programs/2017/MRST_programs/SPE/Results/';
for n=1:N
    fprintf(1,'Main step %d: ',n);
    for m=1:M
        
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
    'dir',dir,'W', W);
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
       x  = explicitTransport(x, G, T/(N*M), rock, fluid, 'wells', W);
       xr  = explicitTransport(xr, G, T/(N*M), rock, fluid, 'wells', W);
       
       pause
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
  
    if n>4, continue, end;
    figure(500)
    % Plot multiphase solution
    subplot(2,4,n);
    contourf(reshape(G.cells.centroids(:,1), G.cartDims),...
        reshape(G.cells.centroids(:,2), G.cartDims), ...
        reshape(x.s(:,1),G.cartDims), [0 cval 1], 'EdgeColor','none');
    hold on;
    
    % Plot corresponding time lines from single-phase solution
    contour(reshape(G.cells.centroids(:,1), G.cartDims),...
        reshape(G.cells.centroids(:,2), G.cartDims), ...
        reshape(tau/T,G.cartDims), a*n/N, '-k','LineWidth',1);
    caxis([0 1]);
    axis equal; axis([0 domain(1) 0 domain(2)]);
    title(sprintf('t=%.2f PVI',n*.2));
    set(gca,'XTick',[],'YTick',[]);
    
    % Plot multiphase solution as function of single-phase time-of-flight
    subplot(2,4,4+n)
    set(gca,'position',get(gca,'position')+[0 .12 0 0]);
    plot(tau(:,1)/tbf,x.s(:,1),'.k','MarkerSize',4);
    set(gca,'XLim',[0 2]); drawnow;
    
    
    
    figure(100)

    if n>4, continue, end;
    
    % Plot multiphase solution
    subplot(2,4,n);
    contourf(reshape(G.cells.centroids(:,1), G.cartDims),...
        reshape(G.cells.centroids(:,2), G.cartDims), ...
        reshape(xr.s(:,1),G.cartDims), [0 cval 1], 'EdgeColor','none');
    hold on;
    
    % Plot corresponding time lines from single-phase solution
    contour(reshape(G.cells.centroids(:,1), G.cartDims),...
        reshape(G.cells.centroids(:,2), G.cartDims), ...
        reshape(tau/T,G.cartDims), a*n/N, '-k','LineWidth',1);
    caxis([0 1]);
    axis equal; axis([0 domain(1) 0 domain(2)]);
    title(sprintf('t=%.2f PVI',n*.2));
    set(gca,'XTick',[],'YTick',[]);
    
    % Plot multiphase solution as function of single-phase time-of-flight
    subplot(2,4,4+n)
    set(gca,'position',get(gca,'position')+[0 .12 0 0]);
    plot(tau(:,1)/tbf,xr.s(:,1),'.k','MarkerSize',4);
    set(gca,'XLim',[0 2]); drawnow;
    
    
    
    
    
    
    
end

%% Plot production curves
% First we plot the saturation in the perforated cell and the corresponding
% water cut, i.e., the fractional flow evaluated in the completion
dt = [0; ones(N*M,1)*T/(N*M)]; t = 1.2*cumsum(dt)/T;
figure;
plot(t,cellfun(@(x) x(2).Sw, wellSols),'--', ...
    t,cellfun(@(x) x(2).wcut, wellSols),'-','LineWidth',1);
legend('Sw in completion','Water cut','Location','NorthWest');
axis([0 max(t) -.05 1.0]);

%%
% Second, we plot the oil rate used in our simulation in units m^3/day
figure;
qOs = cellfun(@(x) abs(x(2).qOs), wellSols);
stairs(t, qOs([2:end end])*day,'LineWidth',1); axis([0 1.2 30 260]);

%%
% Last, we plot the cumulative oil production computed from the well
% solution and compare this with the amount of extracted oil derived from a
% mass-balance computation (initial oil in place minus current oil in
% place). We also include a horizontal line indicating the initial oil in
% place, and a straight line showing the amount of oil we would have
% extracted if oil was produced at the constant initial rate
figure
plot(t,cumsum(bsxfun(@times, abs(cellfun(@(x) x(2).qOs, wellSols)), dt)));
hold on;
plot(t,oip(1)-oip,'-.','LineWidth',3);
plot([0 1.2],oip([1 1]),'-k',t,min(t*oip(1),oip(1)),'-k', ...
    t(W(2).bt+[0 0]),[0 oip(1)],'--k');
hold off; axis tight; axis([0 max(t) 0 1.05*oip(1)]);


%% Plotting with GUI from ad-core
mrstModule add ad-core
plotWellSols(wellSols,cumsum(dt))

%% BS
%% Plot production curves
% First we plot the saturation in the perforated cell and the corresponding
% water cut, i.e., the fractional flow evaluated in the completion
dt = [0; ones(N*M,1)*T/(N*M)]; t = 1.2*cumsum(dt)/T;
figure;
plot(t,cellfun(@(xr) xr(2).Sw, wellSols),'--', ...
    t,cellfun(@(xr) xr(2).wcut, wellSols),'-','LineWidth',1);
legend('Sw in completion','Water cut','Location','NorthWest');
axis([0 max(t) -.05 1.0]);

%%
% Second, we plot the oil rate used in our simulation in units m^3/day
figure;
qOs = cellfun(@(xr) abs(xr(2).qOs), wellSols);
stairs(t, qOs([2:end end])*day,'LineWidth',1); axis([0 1.2 30 260]);

%%
% Last, we plot the cumulative oil production computed from the well
% solution and compare this with the amount of extracted oil derived from a
% mass-balance computation (initial oil in place minus current oil in
% place). We also include a horizontal line indicating the initial oil in
% place, and a straight line showing the amount of oil we would have
% extracted if oil was produced at the constant initial rate
figure
plot(t,cumsum(bsxfun(@times, abs(cellfun(@(xr) xr(2).qOs, wellSols)), dt)));
hold on;
plot(t,oipr(1)-oipr,'-.','LineWidth',3);
plot([0 1.2],oipr([1 1]),'-k',t,min(t*oipr(1),oipr(1)),'-k', ...
    t(Wr(2).bt+[0 0]),[0 oipr(1)],'--k');
hold off; axis tight; axis([0 max(t) 0 1.05*oipr(1)]);


%% Plotting with GUI from ad-core
mrstModule add ad-core
plotWellSols(wellSolsr,cumsum(dt))







