close all
clear all
clc

%% Set up grid and petrophysical data
% We use a Cartesian grid of size nx-by-ny with homogeneous petrophysical
% data: permeability of 100 mD and porosity of 0.2.

%pw=1;
sz=5;
% w1=10*2^pw;
% w2=2*w1;
nxi=1;
nyi=1;
nx=sz;
ny=sz;
iteration=500;
e=1;
num='5';
Lx=sz;
Ly=sz;

for k = 7
    k
tol = 10^-k;

for per=[1]
   
close all
per
%Create the directory
dir='/mnt/sda2/cortes/Results/Solvers/17_09_26/';
%dir='/dev/media/Sphinx/Doctorado_Delft/marzo16/2016_03/Heterogeneous/POD/';
%dir='/dev/media/Sphinx/Doctorado_Delft/Research_programs/2016_1_report/2016_02/heterogeneous/matrixE/';
folder=['size_'   num2str(sz) 'perm_' num2str(per) '_5wells' num2str(k)  ];
mkdir([dir], folder)
dir = [dir folder '/'];
G = cartGrid([sz, sz, 1], [sz, sz, 1]);
G = computeGeometry(G);
% Disable gravity
gravity off

% Set up uniform permeability and constant porosity
rock.perm = ones(G.cells.num, 1)*1;
%inhomogeneus permeability
 lsize=round(sz*sz/sz); 
 
% for i=2:2:6
%  %for i=2
%  if i==2
%      init1=lsize*(i-i/4);
%   rock.perm(1+init1:init1+lsize)  = repmat(10^(-per)*milli*darcy(), [lsize, 1]);
%  else
%      init=init1+lsize*(i-2);
%      rock.perm(init+1:init+lsize)  = repmat(10^(-per)*milli*darcy(), [lsize, 1]);
%  end
%  end
z1=zeros(sz*sz,sz*sz);
nzx=0;
nzy=0;

for i=1:sz

    Z((1:sz)+sz*(i-1),i)=1;
end

%Z=eye(sz*sz);

for i=1:2:sz
 rock.perm(1+lsize*(i-1):lsize*i)  = repmat(10^(-per), [lsize, 1]);
end




z1s=size(Z)


ysz=sz/8;
hy=0;
% %for h=1:sz;
%     for i=1:2:8 
%        
%         for k=1:ysz
%             Z(ysz*(i-1)+1:ysz*i,ysz*(i-1)+1:ysz*i)=1;
%     end
%     end
% %end
% size(Z)
% spy(Z)


 %Visualize the plot
figure
  hp=plotCellData(G, rock.perm);
 
 view(0,90) 
 %return
%  file='perm';
%  axis equal tight
%  B=[dir   file  '.fig'];
%  saveas(hp,B)
%  B=[dir   file   '.jpg'];
%  saveas(hp,B) 
%return
%  break
% rock.perm(1:G.cells.num)  = repmat(10*milli*darcy(), [128, 1]);
rock.poro = ones(G.cells.num, 1)*0.2;
% rock.perm = repmat(100*milli*darcy, [G.cells.num, 1]);
% rock.poro = repmat(0.5            , [G.cells.num, 1]);



%% Compute half transmissibilities
% All we need to know to develop the spatial discretization is the reservoir
% geometry and the petrophysical properties. This means that we can compute
% the half transmissibilities without knowing any details about the fluid
% properties and the boundary conditions and/or sources/sinks that will
% drive the global flow:
hT = simpleComputeTrans(G, rock);

%% Fluid model
% When gravity forces are absent, the only fluid property we need in the
% incompressible, single-phase flow equation is the viscosity. However, the
% flow solver is written for general incompressible flow and requires the
% evaluation of a fluid object that can be expanded to represent more
% advanced fluid models. Here, however, we only use a simple fluid object
% that requires a viscosity and a density (the latter is needed when gravity
% is present)3 bar
gravity reset off
fluid = initSingleFluid('mu' , 1, ...
                        'rho', 1);
display(fluid);
%Boundary conditions
% bc  = pside([], G, 'YMin',3.*barsa());
% bc  = pside(bc, G, 'YMax',0.*barsa());

%% Source terms
% The number of wells can be changes.
% pv  = sum(poreVolume(G,rock));
%% Change number of solution, 1=complete problem, 1,2,... different wells



%% Define well locations, 15 wells

for s=5
    well(1)=2;
    well(2)=-2;
 xi=ones(nx*ny,1)* barsa();

wtype    = {'bhp', 'bhp'};
%wtarget  = [well(1),   well(2)] .* barsa();
wtarget  = [well(1),   well(2)];
wrad     = [0.125, 0.125] .* meter;
wloc     = [  nxi,   nx ;
              nyi,   ny ];
wname    = {'W1', 'W2'};
sgn      = [ 1 ,  -1 ];

W = [];        
for w = 1 : numel(wtype),
   W = verticalWell(W, G, rock, wloc(1,w), wloc(2,w), [], ...
                    'Type', wtype{w}, 'Val', wtarget(w), ...
                    'Radius', wrad(w), 'Name', wname{w}, ...
                    'Sign', sgn(w), 'InnerProduct', 'ip_tpf');
end
nw = length(well);
sz2=size(Z,2);
for nd=Lx*Ly+1:Lx*Ly+nw  
    Z(nd,nd-Lx*Ly)=1;  
end
size(Z)
figure
  hp=plotCellData(G, rock.perm);
hp=plotWell(G, W);
view(0,90)  
%break
sol = initState(G, W, 0);
%load A1;

% Reference TPFA
  mrstModule add agmg
  %mrstModile add PCG_ICSol
   %  solver = AGMGSolverAD('tolerance', 1e-5);
   %solver = GMRES_ILUSolverAD('tolerance', 1e-5);
  cn=0;
%solver = PCG_ICSolverAD_cn('tolerance', tol,'maxIterations', 1000);
solver = PCG_ICSolverAD_cn('tolerance', tol,'maxIterations', 1000,'cn',cn);
%solver = DPCG_ICSolverAD('tolerance', tol,'maxIterations', 1000, 'Z',Z);
% solver = BackslashSolverAD();
% pressureSolver = BackslashSolverAD();
% linsolve = LinearSolverAD('ellipticSolver', pressureSolver);
%linsolver = LinearSolverAD('ellipticSolver', pressureSolver);

linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
tic
psolve = @(state) incompTPFA_g_o(state, G, hT, fluid, 'wells', W,'MatrixOutput',true,'LinSolve', linsolve_p);
% psolve = @(x,p0) ...
%      incompTPFA_g_o(x, G, T, fluid, 'wells', W,'LinSolve', linsolve_p);
[sol,report]= psolve(sol);
toc
p=sol.pressure;
 A=sol.A(1:G.cells.num,1:G.cells.num);
 b=sol.rhs(1:G.cells.num);
 size(A)
% clf;
xb=A\b;
resb=p-xb;
%[xd,flag,res,its]=DICCG_01_25_2(A,b,Z,tol,1000);
%resd=xd-xb;
figure(s+per+200)
[ht]=plotingsolution(G,W,'DICCG', p,1) ;
colorbar
 [ht]=plotingsolution(G,W,'backslash',xb,2);
 colorbar
figure(s+per+300)

plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1);
plotCellData(G,xb-p)
view(0,90)

axis equal tight; colormap(jet(128));
title('bs-DICGCG');
 colorbar
%  figure(s+per+300)
%  [h1]=plotingsolution(G,W,'DICGCG_m',xd-xb,1);
% colorbar




end
end
end
figure; plot(report.residual,'o')
