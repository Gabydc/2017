close all
clear all
clc

%% Set up grid and petrophysical data
% We use a Cartesian grid of size nx-by-ny with homogeneous petrophysical
% data: permeability of 100 mD and porosity of 0.2.
%%If we want to save the files/graphs, it's neccesary to set savefiles==1;
savefiles=1;
%Select a linear solver 1.AGMG, 2.GMRES, 3.PCG_IC, 4.DPCG_IC, other Backslash
lsolver=3;
% Select transport solver 1. Explicit, 2. Implicit
tsolver=2;
% We define the size of the reservoir and the inicial points of the grid
% sz is the number of cells, is the same for the two dimensions
sz=32;
% If we want to compute deflationvectors def==1
%Number of deflation vectors
dv=10;
podv{1}=[1:10];
podv{2}=[6:10];
for i=1:7
f(i) = figure('visible','off');
end    
        % If we use wells, we can define here the position of the wells
        % w1=10*2^pw;
        % w2=2*w1;
        nxi=1;
        nyi=1;
        nx=sz;
        ny=sz;
        %We define the length of the reservoir (physical dimensions)
        Lx=sz;
        Ly=sz;
        
        % We define the maximum of iterations for the liner solver
        maxIterations=500;
        % We define the tolerance of the linear solver
        k=7;
        tol = 10^-k;
        % We can change the permeability one layer is 10*milli*darcy(), the second
        % is 10*milli*darcy()*10^(-per)
        x = linspace(0, 1, 11) .';
        y = linspace(1, 0, 11) .';
        cap_scale = 10;
        pc_form = 'nonwetting';
        [kr, pc]  = tabulatedSatFunc([x, x.^2, y.^2, y.*cap_scale*barsa]);
        props = constantProperties([   1,  10] .* centi*poise, ...
            [1000, 700] .* kilogram/meter^3);
        
        cp=0;
        for def=[0 ]
   
    % If we want to use POD pod==1
    
    for pod = [0 ]
       
       dpod=[];  
       if (def == 0)  && (pod == 1)
           pod=0;
    continue
    end

            for rr=1:2
                dpod=podv{rr};
                if def== 0
                    dpod=[];
                end
 
    
        
          
         
        for per=[1]
            
            close all
            clear Z
            %Create the directory
            dir='/mnt/sda2/cortes/Results/17_03/two_phases/28/';
            
            folder=[ '10-' num2str(k) '_' num2str(sz) 'perm_' num2str(per) 'cp' num2str(cp)];
            mkdir([dir], folder)
            dir1 = [dir folder '/'];
            
            folder=[  'def_' num2str(def) '_pod_' num2str(numel(dpod)) ];
            mkdir([dir1], folder)
            dir2 = [dir1 folder '/'];
            %handler2 = ResultHandler('writeToDisk', true,'dataDirectory',dir,'dataFolder','results');
            %Create the grid
            domain = [250 250 1];
            G = cartGrid([sz, sz, 1], domain);
            G = computeGeometry(G);
            % Disable gravity
            gravity off
            
            % Set up uniform permeability and constant porosity
            rock.perm = ones(G.cells.num, 1)*1*milli*darcy();
            rock.poro = ones(G.cells.num, 1)*0.2;
            % Create layers of permeability
            %lsize=round(sz*sz/sz);
            %     for i=1:2:sz
            %         rock.perm(1+lsize*(i-1):lsize*i)  = repmat(10^(-per), [lsize, 1])*10*milli*darcy();
            %     end
            lsz=round(sz/8);
            
            for i= 1:2:8
                rock.perm(1+(i-1)*lsz*nx:i*lsz*nx)  = 1*10^(-per)*milli*darcy();
            end
            f(1)=figure(1);
            plotCellData(G, rock.perm);
            view(0,90), colormap(jet), axis equal tight
  
            
            % Create the subdomain-based deflation vectors
            %     for i=1:sz
            %         Z((1:sz)+sz*(i-1),i)=1;
            %     end
            %     hz=figure;
            %     spy(Z);
            %
            %     if savefiles==1
            %         file='Z';
            %         axis tight
            %         B=[dir   file  '.fig'];
            %         saveas(hz,B)
            %         B=[dir   file   '.jpg'];
            %         saveas(hz,B)
            %         %Z=eye(sz*sz);
            %     end
            
            
            %% Compute half transmissibilities
            % All we need to know to develop the spatial discretization is the reservoir
            % geometry and the petrophysical properties. This means that we can compute
            % the half transmissibilities without knowing any details about the fluid
            % properties and the boundary conditions and/or sources/sinks that will
            % drive the global flow:
            hT = simpleComputeTrans(G, rock);
            
            %% Fluid model
            
            % We set up a two-phase fluid. Viscosity, density is set for oil and water
            
            gravity off
            verbose = false;
            if cp==0
                fluid = initSimpleFluid('mu' , [   1,    10] .* centi*poise     , ...
                    'rho', [1000, 700] .* kilogram/meter^3, ...
                    'n'  , [   2,    2]);
            else
                
                fluid = struct('properties', props                  , ...
                    'saturation', @(x, varargin)    x.s  , ...
                    'relperm'   , kr                     , ...
                    'pc'        , @(x, varargin) pc(x.s));
                xDummy   = initState(G, [], [0, 1]);
                xDummy.s = linspace(0, 1, numel(xDummy.s))'; ...
                    pc = convertTo(fluid.pc(xDummy), barsa);
                
                hcp=figure(232);
                plot(xDummy.s, pc);
                xlabel('s_w'); ylabel('pc [bar]');
                title('Capillary pressure curve')
                
            end
            
            s = linspace(0, 1, 1001).'; kr = fluid.relperm(s);
            f(2) = figure(2);
            plot(s, kr); legend('kr_1', 'kr_2'),title(['rel perm']), axis equal tight
                   
 
            
            %Boundary conditions. We set boundary conditions, waterflooding
            
            pv = poreVolume(G, rock);
            %injRate = -sum(pv)/(500*day);
            injRate = -0.1*meter^3/day;
            %     bc = fluxside([], G, 'ymin', -injRate, 'sat', [1, 0]);
            %     bc = pside(bc, G, 'ymax', 0*barsa, 'sat', [0, 1]);
            W = addWell([],G, rock, 1, 'Type', 'bhp', ...
                'Val', 100*barsa, 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
            W = addWell(W,G, rock, G.cells.num, 'Type', 'bhp', ...
                'Val', 0, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);
            
            %% Set pressure solver
            
            
            switch lsolver
                case 1
                    mrstModule add agmg
                    solver = AGMGSolverAD('tolerance', tol);
                case 2
                    solver = GMRES_ILUSolverAD('tolerance', tol);
                case 3
                    solver = PCG_ICSolverAD('tolerance', tol,'maxIterations', maxIterations);
                case 4
                    solver = DPCG_ICSolverAD('tolerance', tol,'maxIterations', maxIterations, 'Z',Z);
                otherwise
                    solver = BackslashSolverAD();
            end
            %S  = computeMimeticIP(G, rock, 'Verbose', verbose,'InnerProduct','ip_tpf');
            fn = @(A, b) solver.solveLinearSystem(A, b);
            psolve = @(state) incompTPFA_g_o(state, G, hT, fluid, 'wells', W, 'MatrixOutput',true,'LinSolve', fn);
            %% Define transport solver
            
            switch tsolver
                case 1
                    % Explicit tranport solver
                    tsolve  = @(state, dT, fluid) explicitTransport(state, G, dT, rock, fluid, 'wells', W, 'verbose', false);
                case 2
                    % Implicit transport solver: try with one time step
                    tsolve  = @(state, dT, fluid) implicitTransport(state, G, dT, rock, fluid, 'wells', W, 'Verbose', false);
            end
            
            %% Initiate pressure solver
            rSol = initState(G, [], 100*barsa, [0 1]);
            
            heading = [' saturation'];
            %heading = [num2str(t),  ' days'];
            f(3) = figure(3);
            subplot('position',[0.001, 0, 0.5, 0.9]); cla
            plotCellData(G, rSol.s(:,1));
            caxis([0 1]), view(0,90), axis equal off, title(['Water' heading]), colormap(jet),colorbar
            subplot('position',[0.501, 0, 0.5, 0.9]), cla
            plotCellData(G, rSol.s(:,2));
            caxis([0 1]), view(0,90), axis equal off, title(['Oil' heading]), colormap(jet), colorbar

            
            
    
            
  
            
            
            
            %% Solve initial pressure in reservoir
            rSol = psolve(rSol);
        
            heading = [' pressure'];
            %heading = [num2str(t),  ' days'];
            f(4) = figure(4);
            subplot('position',[0.001, 0, 0.5, 0.9]); cla
            plotCellData(G, rSol.s(:,1));
            caxis([0 1]), view(0,90), axis equal off, title(['Water' heading]), colormap(jet),colorbar
            subplot('position',[0.501, 0, 0.5, 0.9]); cla
            plotCellData(G, rSol.s(:,2));
            caxis([0 1]), view(0,90), axis equal off, title(['Water' ' flux']), colormap(jet), colorbar


            
            
            
            %% Transport loop
            T      = 900*day();
            %T      = 100*day();
            
            dT     = T/180;
            dTplot = 300*day();  % plot only every 100th day
            N      = fix(T/dTplot);
            pv     = poreVolume(G,rock);
            t      = 0;
            plotNo = 1;
            H1 = 'Water Sat - '; H2 = 'Oil Sat - '; H3 = 'W P - ';
            e = []; p_org = []; p_pc = [];
            %figure;
            %dT=20;
            ts=0;
            podi=0;
            while t < T,
                
                ts=ts+1;
                % Increase time and continue if we do not want to plot saturations
 
                if t==0
                    p0=rSol.pressure;
                    cmin=1*min(p0);
                    cmax=3*max(p0);
                    cm=9.3252e+06;
                end

                % TRANSPORT SOLVE
                [rSol,treport(ts)]   = tsolve(rSol, dT, fluid);
                %rSol_pc = tsolve(rSol_pc, dT, fluid_pc);
                
                % Check for inconsistent saturations
                %s = [rSol.s(:,1); rSol_pc.s(:,1)];
                s = [rSol.s(:,1)];
                assert(max(s) < 1+eps && min(s) > -eps);
                
                
         
                
                
                if  def==0
                    [rSol,preport(ts)]    = psolve(rSol);
                else
                    Zp(:,ts)=rSol.pressure;
                    if ts <dv+1
                        % Update solution of pressure equation. If backslash is used, there
                        % is no report
                        [rSol,preport(ts)]    = psolve(rSol);
                        %rSol_pc = psolve(rSol_pc, fluid_pc);
                        
                    else
                        
                        podi=podi+1;
                        Z=Zp(:,podi:podi+dv-1);
                        
                        % pause
                        if pod==1
                            
                            [U,S]=defpodf_Dt(Z,dir2,dv,ts,dTplot/day());
                            Z=U(:,dpod);
                        end
                        
                        for i=1:numel(W)
                            
                            Z(nx*ny+i,:)=0;
                        end
                        
                        solver = DPCG_ICSolverAD('tolerance', tol,'maxIterations', maxIterations, 'Z',Z);
                        fn = @(A, b) solver.solveLinearSystem(A, b);
                        psolve = @(state) incompTPFA_g_o(state, G, hT, fluid, 'wells', W, 'MatrixOutput',true,'LinSolve', fn);
                        [rSol,preport(ts)]    = psolve(rSol);
                    end
                end
                
%                 hst=figure(13);
%                 heading = [' saturation'];
%                 %heading = [num2str(t),  ' days'];
%                 
%                 subplot('position',[0.001, 0, 0.5, 0.9]), cla
%                 plotCellData(G, rSol.s(:,1));
%                 caxis([0 1]), view(0,90), axis equal off, title(['Water' heading]),colorbar
%                 subplot('position',[0.501, 0, 0.5, 0.9]), cla
%                 plotCellData(G, rSol.pressure);
%                 view(0,90), axis equal off, title(['Pressure']), colormap(jet), colorbar
%                 %pause
                
                
                
                
                
                t = t + dT;
                
                if ( t < plotNo*dTplot && t<T) continue, end
                % Plot saturation
                heading = [num2str(convertTo(t,day)),  ' days'];
                %heading = [num2str(t),  ' days'];
                r = 0.01;
                f(5) = figure(5);
                subplot('position',[(plotNo-1)/N+r, 0.50, 1/N-2*r, 0.4]); cla
                plotCellData(G, rSol.s(:,1));
                caxis([0 1]), view(0,90), axis equal off, title([H1 heading]), colormap(jet)
                subplot('position',[(plotNo-1)/N+r, 0.02, 1/N-2*r, 0.4]); cla
                plotCellData(G, rSol.s(:,2));
                caxis([0 1]), view(0,90), axis equal off, title([H2 heading]), colormap(jet)
                            
                
                %heading = [num2str(t),  ' days'];
                r = 0.01;
                
                if t==dTplot;
                    heading = [num2str(convertTo(dT,day)),  ' days'];
                    f(6) = figure(6);
                    subplot('position',[(plotNo-1)/(N+1)+r, 0.1, 1/(N+1)-r, 0.9]), cla
                    plotCellData(G, p0);
                    view(0,90), %colorbar,
                    axis equal tight off, colormap(jet), title([H3 heading]),caxis([cmin cm])
                end
                %pause
       
               
                heading = [num2str(convertTo(t,day)),  ' days'];
                p=rSol.pressure;
                f(6) = figure(6);
                subplot('position',[(plotNo)/(N+1)+r, 0.1, 1/(N+1)-r, 0.9]); cla
                plotCellData(G, p);
                view(0,90), %colorbar,
                axis equal tight off, colormap(jet), title([H3 heading]),caxis([cmin cm])
                plotNo = plotNo+1;
            end
            
            
            
            p=rSol.pressure;
            
            A=rSol.A(1:G.cells.num,1:G.cells.num);
            b=rSol.rhs(1:G.cells.num);
            
            % clf;
            xb=A\b;
            resb=p-xb;
            %[xd,flag,res,its]=DICCG_01_25_2(A,b,Z,tol,1000);
            %resd=xd-xb;
            
            figure(per+1200)
            title('Pressure')
            [ht]=plotingsolution_bc(G,'DICCG', p,1) ;
            colorbar
            [ht]=plotingsolution_bc(G,'backslash',xb,2);
            colorbar
            
            figure(per+1300)
            
            plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1);
            plotCellData(G,xb-p)
            view(0,90)
            
            axis equal tight; colormap(jet(128));
            title('bs-DICGCG');
            colorbar
            
            for j=1:ts
                its(j)=preport(j).iterations;
            end
            f(7) = figure(7);
            if def==1
               
                plot(1:dv,its(1:dv),'r*');
                hold on
                plot(dv+1:ts,its(dv+1:ts),'bp');
                hold on
                legend('ICCG','DICCG');
                axis square
            else
                plot(1:ts,its(1:ts),'r*');
                legend('ICCG');
                axis square
            end
            ylabel('Number of iterations','FontSize',16)
            xlabel('Time step ','FontSize',16)
            if def==0
                title(['Iterations'],'FontSize',16)
            else if (def==1)&&(pod==0)
                    title(['Iterations, ' num2str(dv) ' deflation vectors'],'FontSize',16)
                else
                    title(['Iterations, ' num2str(dv) ' snapshots, ' num2str(numel(dpod)) ' POD vectors '],'FontSize',16)
                end
            end
            if savefiles==1
                
                file{1} = 'Permeability';
                file{2} = 'RelPerm';
                file{3} = 'ISat';
                file{4} = 'IPressure';
                file{5} = 'Sat';
                file{6} = 'Pressure';
                file{7} = 'Iterations';
                if cp==1
                    
                    filecp='cp';
                    axis equal tight
                    B=[dir2   filecp  '.fig'];
                    saveas(hcp,B)
                    B=[dir2   filecp   '.jpg'];
                    saveas(hcp,B)
                end
                for i=1:7
                   
                    axis equal tight
                    axis square
                    %i
                B=[dir2   file{i}  '.fig'];
                savefig(f(i),B)
                 B=[dir2   file{i}   '.jpg'];
                 saveas(f(i),B)
                end
                 filepr=['preport'];
                filename=[dir2 filepr];
                save(filename,filepr)
                filepr=['treport'];
                filename=[dir2 filepr];
                save(filename,filepr)
                filetx = ['results.txt'];
                saveres(dir1,filetx,def,pod,dpod,per,ts,dv,preport)
                addplots(dir1,dir2,file,'plots.txt')
              end  
   
            
        end
    end
        end
        end
  
        
figure(6)





