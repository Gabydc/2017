%% Single phase flow simulation using AD
% This example goes through the steps of setting up a single-phase
% simulation with a single horizontal well using the automatic
% differentiation framework.
close all
clear all
clc
% Required modules
try
    require ad-fi
catch
    mrstModule add ad-fi
end
%save graphs
ss=1;
%solve with deflation if not
def=0;
defvect=5;
% Setup 10x10x10 grid of 200x200x50 m model.
nx = 35;    ny = 35;    nz = 1;
Dx = 35;   Dy = 35;   Dz = 50;
G = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
%G = cartGrid([nx, ny, nz]);
G = computeGeometry(G);

%% Setup rock properties
%

% Assume homogeneous/isotropic rock.
permX = 3*milli*darcy;
poro  = 0.3;
rock.perm = repmat(permX, [G.cells.num, 1]);
rock.poro = repmat(poro , [G.cells.num, 1]);
sz=nx;
lsz=round(sz/5);
per=0;
for i= 1:2:5
    rock.perm(1+(i-1)*lsz*nx:i*lsz*nx)  = 30*10^(-per)*milli*darcy();
end
dir='/home/wagm/cortes/Localdisk/Results/16_06/17/';
%dir='/dev/media/Sphinx/Doctorado_Delft/Results/16_06/06/';
folder=['size_'   num2str(nx) 'perm_' num2str(per) '_5wells_' num2str(defvect) '_defvect' num2str(def) 'c01'];
mkdir([dir], folder)
dir = [dir folder '/'];


 hp=plotCellData(G, rock.perm);

% Set rock compressibility:
cr = 1e-6/barsa;
c    = 1e-3/barsa;

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
p_r  = 200*barsa;

% Finally, the pressure-dependent function for pore-volumes becomes:

pv   = @(p) pv_r .* exp( cr .* (p - p_r) );

%% Fluid (oil) properties
% We Assume constant viscosity:
mu   = 5*centi*poise;
% Assume that the oil compressibility can be approximated as constant in
% the reservoir:

% When the compressibility is constant, the fluid density becomes a
% function of pressure given by the differential equation
%                 c = (d rho/d p)/rho
% which results in
%                 rho(p) = rho_r e^( c(p-p_r) )
% in which rho_r is the reference fluid density at reference pressure p_r.
% We assume that rho_r = 850 kg/m^3 at p_r = 200 Bar:

rho_r = 850*kilogram/meter^3;

% Finally define the pressure dependent function for rho:
rho   = @(p) rho_r .* exp( c .* (p - p_r) );

% Computing surface volume rates requires fluid density at surface
% conditions. We assume it is 750 kg/m^3:
rhoS = 750*kilogram/meter^3;

%% Single horizontal well in J direction
% We consider a well with 8 connections in the J direction.

wloc=[1  nx    1  nx   ceil(nx/2);
      1   1   ny  ny   ceil(ny/2);
      1   1    1   1            1];
  
W = [];
for i=1:5
I = wloc(1,i);
J = wloc(2,i);
K = wloc(3,i);
cellInx = sub2ind(G.cartDims, I, J, K);
W = addWell(W, G, rock, cellInx, 'Name', 'producer', 'InnerProduct', 'ip_tpf', 'sign', -1);
end

% Plotting
figure(520)
hp=plotCellData(G, rock.perm);
plotWell(G, W);
axis off; 
color 'w';
camproj perspective;
view(0,90)



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
%ad_values is a vector containing the resulting initialized adi variables
%[output{1}, output{2}, output{3},output{4}, output{5}, output{6},output{7}, output{8},output{9},output{10},output{11}]=
%[     p_ad,    bhp_ad,   bhp_ad2,  bhp_ad3,   bhp_ad4,   bhp_ad5,    qS_ad,    qS_ad2,   qS_ad3,    qS_ad4,    qS_ad5]
[output{1}, output{2}, output{3},output{4}, output{5}, output{6},output{7}, output{8},output{9},output{10},output{11}]=...
    initVariablesADI(inputv{1},inputv{2},inputv{3},inputv{4},inputv{5},inputv{6},inputv{7},inputv{8},inputv{9},inputv{10},inputv{11});
%return
%[p_ad, bhp_ad, bhp_ad2, qS_ad, qS_ad2] = initVariablesADI(p_init, p_init(wc(1)), p_init(wc2(1)), 0, 0);


%[p_ad,bhp_ad1,bhp_ad2,qS_ad1,qS_ad2]=initVariablesADI(p_init,p_init(wc(1)),p_init(wc(tb)),0,0);
% For convenience, create indices to variables when stacked:
%Ix(1) = pIx, Ix(2)=bhpIx,  Ix(3)=qSIx,Ix(4)=bhpIx2 Ix(5)=qSIx2,
%Ix(6)=bhpIx3, Ix(7)=qSIx3 Ix(8)=bhpIx4, Ix(9)=qSIx4, Ix(10)=bhpIx5,
%Ix(11)= qSIx5 ;
% pIx = 1:G.cells.num;
% bhpIx  = G.cells.num + 1; qSIx  = G.cells.num + 6;
% bhpIx2 = G.cells.num + 2; qSIx2 = G.cells.num + 7;
% bhpIx3 = G.cells.num + 3; qSIx3 = G.cells.num + 8;
% bhpIx4 = G.cells.num + 4; qSIx4 = G.cells.num + 9;
% bhpIx5 = G.cells.num + 5; qSIx5 = G.cells.num + 10;
Ix{1} = 1:G.cells.num;
for i=1:10
    Ix{i+1}=G.cells.num+i;
end

%% Set up simulation parameters
%

numSteps = 52;
totTime  = 156*day;
dt       = totTime / numSteps;
% Tolerance and maximum number of iterations for the Newton solver.
tolNR      = 1e-7;
maxits   = 10;
miter=zeros(numSteps,maxits);
% Save output in array 'sol'
sol = repmat(struct('time', [], 'pressure', [], 'bhp', [], 'bhp2', [], ...
    'qS', [], 'qS2', []), ...
    [numSteps + 1, 1]);

sol(1).time      = 0;
sol(1).pressure  = double(output{1});
sol(1).bhp       = double(output{2});
sol(1).bhp2      = double(output{3});
sol(1).bhp3      = double(output{4});
sol(1).bhp4      = double(output{5});
sol(1).bhp5      = double(output{6});
sol(1).qS        = double(output{7});
sol(1).qS2       = double(output{8});
sol(1).qS3       = double(output{9});
sol(1).qS4       = double(output{10});
sol(1).qS5       = double(output{11});
% Set up plot
 f = [figure(1), figure(2)];

clf
set(f(2), 'Color', 'w')
subplot(2,2,1),
plotCellData(G, convertTo(p_init, barsa));
title('pressure [bar]', 'EdgeColor', 'w');
%colorbar;
view(0,90);
camproj perspective

subplot(2,2,[3,4]);
xlim([0, convertTo(totTime,day)]);
title('Surface volume rate [m^3/day]');
hold on


%% Main simulation
% We solve the equations implicitely. At each time step, the equations are assembled and
% the automatic differentiation framework takes care automatically of the computation of
% the Jacobian.

t = 0;
step = 0;
counter=0;
%initialize counter for the snapshots
sv=0;
dv1=0;
nDigits = floor(log10(maxits)) + 1;
text = [dir 'iterations.txt'];
fileID = fopen(text,'a');
fprintf(fileID,'Time step & Newton iteration & # Iterations solver \\\\ \n');
fclose(fileID);
% fprintf(fileID,' %6s', '%6s ', '%6.2e ', step, nit, iter ' \\');
upd=sol(1).pressure;


while t < totTime,
    pf(1)=sol(step+1).pressure(1);
    nns=1;
    for i=1:nx-1
        nns=nns+1;
        ns=1+nx*i+i;
        pf(nns)=sol(step+1).pressure(ns);
    end
    set(0, 'CurrentFigure', f(2));
    subplot(2,2,2),  %caxis([120, 205])
    plot(convertTo(pf, barsa));
    xlim([1, nx]);
    title('Pressure diagonal [bar]');
    hold on
    subplot(2,2,[3,4]);
    xlim([0, convertTo(totTime,day)]);
    title('Surface volume rate [m^3/day]');
    hold on
    t = t + dt;
    
    step = step + 1;
    fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
        step, convertTo(t - dt, day), convertTo(t, day));
  textnr = [dir 'NRiterations.txt'];
  fileID = fopen(textnr,'a');
  fprintf(fileID,'\nTime step %d: Time %.2f -> %.2f days \\\\ \n', ...
        step, convertTo(t - dt, day), convertTo(t, day));
  fclose(fileID);
    
    % Newton loop
    resNorm = 1e99;
    p0  = double(output{1}); % Previous step pressure
    nit = 0;
    
    if mod(step,5)==0 && step>5;
        counter=counter+1;
    end
    leginf= cell(10,0);
    while (resNorm > tolNR) && (nit < maxits)
        % Create equations:
        %eqs = cell([3, 1]);
        
        eqs{1} = pressureEq(output{1}, p0, dt);
        % Add well contributions in perforated cells:
        eqs{1}(W(1).cells) = eqs{1}(W(1).cells) - wellRates{1}(output{1}, output{2});
        eqs{1}(W(2).cells) = eqs{1}(W(2).cells) - wellRates{2}(output{1}, output{3});
        eqs{1}(W(3).cells) = eqs{1}(W(3).cells) - wellRates{3}(output{1}, output{4});
        eqs{1}(W(4).cells) = eqs{1}(W(4).cells) - wellRates{4}(output{1}, output{5});
        eqs{1}(W(5).cells) = eqs{1}(W(5).cells) - wellRates{5}(output{1}, output{6});
        
        % Final equation is prescribed bhp at producers (injector)
        eqs{2} = output{2} - 100*barsa;
        eqs{3} = output{3} - 100*barsa;
        eqs{4} = output{4} - 100*barsa;
        eqs{5} = output{5} - 100*barsa;
        eqs{6} = output{6} - 600*barsa;
        % Sum of wellrates should equal total rate:
        eqs{7} = output{7} - sum(wellRates{1}(output{1}, output{2}))/rhoS;
        eqs{8} = output{8} - sum(wellRates{2}(output{1}, output{3}))/rhoS;
        eqs{9} = output{9} - sum(wellRates{3}(output{1}, output{4}))/rhoS;
        eqs{10} = output{10} - sum(wellRates{4}(output{1}, output{5}))/rhoS;
        eqs{11} = output{11} - sum(wellRates{5}(output{1}, output{6}))/rhoS;
        % Concatenate equations and solve:
        eq  = cat(eqs{:});
        J   = eq.jac{1};  % Jacobian
       
        res = eq.val;     % residual
        %resNorm = norm(res);
       

        % upd = -(J \ res); % Newton update
        l=ichol(J);
        e=0;
        if step==1
            xi=ones(length(res),1);
        else
            xi=upd;
        end
        ne=0;
  
        if resNorm > tolNR
        if  def==0
             clear ee
             [upd,ee]=ICCG_t1x(-J,res,xi,200,10^-7,l,dir,tolNR);
            ne=length(ee);
            miter(step,nit+1)=ne;
            
               fileID = fopen(text,'a');
                %fprintf(fileID,'Time step', 'Newton iteration', '# Iterations solver \\');
                fprintf(fileID,' %6d  &   %6d  &  %6d \\\\ \n', step, nit+1, ne);
                fclose(fileID);
        else
            if step <defvect+1
                if nit==0
                    %counter for the snapshots
                    sv=sv+1;
                end
                [upd,ee]=ICCG_t1x(-J,res,xi,200,10^-7,l,dir,tolNR);
                %vector containing the snapshots
                zd(:,sv)=upd;
                %number of iterations of the ICCG method
                ne=length(ee);
                %the number of iterations of ICCG method is stored in a matrix for each timestep 
                miter(step,nit+1)=ne;
                iter=1:ne;
%                            figure(800)
%                             if defvect<defvect+1
%                                 subplot(2,3,step)
%                             else if defvect<11
%                                     subplot(2,5,step)
%                 
%                                 else
%                                     subplot(5,5,step)
%                                 end
%                             end
%                             colors = [1 0 0;0 1 0; 0 1 1; 1 0 1; 0 0 1;1 0 1;0.2 0.5 0.5; 0.5 0.2 1; 1 0.2 0.5; 0.2 1 0.6;];
%                             hc=semilogy(iter,ee,'o','Color',colors(nit+1,:));
%                 
%                             hold on
%                             title( ['Timestep ' num2str(step )] )
%                             leginf{nit+1} = ['N-R ' num2str(nit+1)];
%                            legend(leginf{1:nit+1})
                fileID = fopen(text,'a');
                %fprintf(fileID,'Time step', 'Newton iteration', '# Iterations solver \\');
                fprintf(fileID,' %6d  &   %6d  &  %6d \\\\ \n', step, nit+1, ne);
                fclose(fileID);
                
                
            else
                if step==defvect+1
                    [U,S,hpod1]=defpodf_D(zd,dir);
                end
                z=S(:,1:5);
                [upd,ee]=DICCG_t(-J,res,xi,200,10^-7,z,l);
            ne=length(ee);
            miter(step,nit+1)=ne;
            iter=1:ne;
            end
             end
       
        
        
            file=['miter'];
            filename=[dir file];
            save(filename,file)
            
            
            
            %             if mod(step,5)==0
            %                 figure(500)
            %                 if def==1
            %                     title('Iterations DICCG');
            %                 else
            %                     title('Iterations ICCG')
            %                 end
            %                 subplot(2,5,counter)
            %                 colors = [1 0 0;0 1 0; 0 1 1; 1 0 1; 0 0 1;1 0 1;0.2 0.5 0.5; 0.5 0.2 1; 1 0.2 0.5; 0.2 1 0.6;];
            %                 hd=semilogy(iter,ee,'o','Color',colors(nit+1,:));
            %
            %                 hold on
            %                 title( ['Timestep ' num2str(step )] )
            %                 leginf{nit+1} = ['N-R ' num2str(nit+1)];
            %                 legend(leginf{1:nit+1})
            %             end
%             fileID = fopen(text,'a');
%             %fprintf(fileID,'Time step', 'Newton iteration', '# Iterations solver \\');
%             ne=length(ee);
%             
%             fprintf(fileID,' %6d  &   %6d   & %6d \\\\  \n ', step, nit+1, ne);
%             fclose(fileID);
            % [upd,hc0,hc01,hc02]=ICCG_t(-J,res,xi,iter,10^-7,l,e,step);
        
      
        %legend(leginf{1:nit+1})
        
        
        % Update variables
         nit     = nit + 1;
         resNorm = norm(res);
        for i=1:11
        output{i}.val   = output{i}.val   + upd(Ix{i});
        end
        else
             for i=1:11
        output{i}.val   = output{i}.val   ;
             end
          nit     = nit + 1;
         resNorm = norm(res);
        end
        fprintf('  Iteration %*d:  Res = %.4e\n', nDigits, nit, resNorm);
          fileID = fopen(textnr,'a');
       fprintf(fileID,'  Iteration %*d:  Res = %.4e \\\\ \n', nDigits, nit, resNorm);
       fclose(fileID);
       
    end
    
    if nit > maxits,
        error('Newton solves did not converge')
    else
        sol(step+1).time     = t;
        sol(step+1).pressure = double(output{1});
        sol(step+1).bhp      = double(output{2});
        sol(step+1).qS       = double(output{7});
        sol(step+1).bhp2      = double(output{3});
        sol(step+1).qS2       = double(output{8});
        sol(step+1).bhp3      = double(output{4});
        sol(step+1).qS3       = double(output{9});
        sol(step+1).bhp4      = double(output{5});
        sol(step+1).qS4       = double(output{10});
        sol(step+1).bhp5      = double(output{6});
        sol(step+1).qS5       = double(output{11});
        % Plot evolution
        set(0, 'CurrentFigure', f(2));
        subplot(2,2,1), cla, %caxis([120, 205])
        plotCellData(G, convertTo(sol(step+1).pressure, barsa), 'EdgeColor', 'w');
        
        subplot(2,2,[3,4])
        xaxis = convertTo(sol(step+1).time, day);
        y1 = convertTo(-sol(step+1).qS , meter^3/day);
        y2 = convertTo(-sol(step+1).qS2 , meter^3/day);
        y3 = convertTo(-sol(step+1).qS3 , meter^3/day);
        y4 = convertTo(-sol(step+1).qS4 , meter^3/day);
        y5 = convertTo(-sol(step+1).qS5 , meter^3/day);
        
        hs= plot(xaxis, y1,'x', xaxis,y2,'r*', xaxis,y3,'bo', xaxis,y4,'gp', xaxis,y5,'ms');
        ylabel('Rate [meter^3/day]','FontSize',14)
        xlabel('Time [days]','FontSize',14)
        legend('P1','P2','P3','P4','I1')
        
        drawnow
        
    end
    
    
end


figure(400)
for i=1:4
    
    if def==1
        subplot(2,2,i)
        hn=plot(1:defvect,miter(1:defvect,i),'r*');
        hold on
        subplot(2,2,i)
        hn=plot(defvect+1:numSteps,miter(defvect+1:numSteps,i),'bp');
        hold on
        legend('ICCG','DICCG');
    else
        subplot(2,2,i)
        hn=plot(miter(:,i),'r*');
        legend('ICCG');
    end
    ylabel('Number of iterations','FontSize',14)
    xlabel('Time step ','FontSize',14)
    title(['Iterations N-R' num2str(i) ],'FontSize',14)
end
if ss==1
    file='permeability_field';
    B=[dir  file  '.fig'];
    saveas(hp,B)
    B=[dir  file  '.jpg'];
    saveas(hp,B)
    file='solution';
    B=[dir  file  '.fig'];
    saveas(hs(1),B)
    B=[dir  file  '.jpg'];
    saveas(hs(1),B)
    
    file='iterations_4NR';
    B=[dir  file  '.fig'];
    saveas(hn,B)
    B=[dir  file  '.jpg'];
    saveas(hn,B)
    if def==1
        file='iter_ICCG(eigenvalues)';
        B=[dir  file  '.fig'];
        saveas(hpod1,B)
        B=[dir  file  '.jpg'];
        saveas(hpod1,B)
    end
end