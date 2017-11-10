close all
clear all
clc

%% Set up grid and petrophysical data
% We use a Cartesian grid of size nx-by-ny with homogeneous petrophysical
% data: permeability of 100 mD and porosity of 0.2.

pw=1;
size=32*2^pw;
w1=10*2^pw;
w2=2*w1;
nx=size;
ny=size;
iteration=500;
k=10;
tol=10^-k;
e=0;
for perm=1:2:7

%Create the directory
dir='/home/wagm/cortes/Localdisk/Research_programs/2016_1_report/1_18_results/';
folder=['heterogeneous_' num2str(size)  '_' num2str(perm)];
mkdir([dir], folder)
dir = [dir folder '/'];
G = cartGrid([size, size, 1]);
G = computeGeometry(G);
% Disable gravity
gravity off

% Set up uniform permeability and constant porosity
rock.perm = ones(G.cells.num, 1)*1*milli*darcy;
% %inhomogeneus permeability
lsize=round(size*size/8);
for i=1:2:8
 rock.perm(1+lsize*(i-1):lsize*i)  = repmat(10^(-perm)*milli*darcy(), [lsize, 1]);
end
% % Visualize the plot
%  plotCellData(G, rock.perm)
% 
%  view(0,90) 
% break
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
fluid = initSingleFluid('mu' , 1*centi*poise, ...
                        'rho', 1014*kilogram/meter^3);
display(fluid)
%Boundary conditions
bc  = pside([], G, 'YMin',3.*barsa());
bc  = pside(bc, G, 'YMax',0.*barsa());

%% Source terms
% The number of wells can be changes.
% pv  = sum(poreVolume(G,rock));
%% Change number of solution, 1=complete problem, 1,2,... different wells



%% Define well locations, 4 wells

for s=1:5
    %pause
    well(1:4)=0;
num=num2str(s);
    if s==1
    xi=rand(nx*ny,1);
         file=[dir 'xi'];
        save(file,'xi')
        well(1)=-5;
    else if s==2
        well(2)=5;
    else if s==3
        well(3)=5;
        else if s==4
               well(4)=-5;
            else
%                 xi=rand(nx*ny,1);
%          file=[dir 'xi'];
%         save(file,'xi')
                well(1)=-5;
                well(2)=5;
                well(3)=5;
                well(4)=-5;
            end
        end
        end
    end
%        well
%        pause
       
wtype    = {'bhp', 'bhp', 'bhp', 'bhp'};
wtarget  = [well(1),   well(2),   well(3),   well(4)] .* barsa();
wrad     = [0.125, 0.125, 0.125, 0.125] .* meter;
wloc     = [  w1,   w1,     w2,   w2;
              w1,   w2,     w1,   w2];
wname    = {'W1', 'W2', 'W3', 'W4'};
sgn      = [ 1 ,  1 ,  1 ,  1  ];
W = [];        
for w = 1 : numel(wtype),
   W = verticalWell(W, G, rock, wloc(1,w), wloc(2,w), [], ...
                    'Type', wtype{w}, 'Val', wtarget(w), ...
                    'Radius', wrad(w), 'Name', wname{w}, ...
                    'Sign', sgn(w), 'InnerProduct', 'ip_tpf');
end

% figure
% plotCellData(G, rock.perm);
% plotWell(G, W);
% break
%   view(0,90) 
% 
% file='perm';
%  B=[dir   file def '.fig'];
%  saveas(hj,B)
%  file='perm';
%  B=[dir   file def '.jpg'];
%  saveas(hj,B)
%  break

% Create a initialized state and set initial saturation to phase 1.
sol = initState(G, W, 0);
% 
% % Find transmissibility.
% T = computeTrans(G, rock);

% Reference TPFA
psolve = @(state) incompTPFA(state, G, hT, fluid, 'wells', W,'MatrixOutput',true,'bc',bc);

sol= psolve(sol);
p=sol.pressure;
 A=sol.A(1:G.cells.num,1:G.cells.num);
 b=sol.rhs(1:G.cells.num);
% clf;
figure
plotingsolution(G,W,'s', sol.pressure,1) 
plotWell(G, W);
view(30,50)  
%%Sol ICCG
if s==1
l=ichol(A);
 file=[dir 'l'];
 save(file,'l')
else
    file=[dir 'l'];
load(file)
end
fc=200;
nn=0;
def=0;
def=['tol-' num2str(k) '_' num2str(s) 'snapshots'];
text = [dir 'iter' def '.txt'];
text1 = [dir 'cond' def '.txt'];

%[x11]=ICCG_1(A,b,xi,iteration,tol,l); 
[x11,iter2,e2,hl2,h2,h21]=ICCG_1(A,b,xi,iteration,tol,l,e,text1,nn,dir,def,fc);
 x(:,s)=x11;
 figure
 [h1]=plotingsolution(G,W,'CG',x(:,s),2);

% colorbar
  
% [i j k] = ind2sub(G.cartDims, 1:G.cells.num);
% clf;
% % Plot the grid
% plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1)
% plotCellData(G, )
% colorbar
% Plot the wells
% figure

% view(30,50)
%  figure
% plotingsolution(G,W,'s',sol.pressure,1);

 file=[dir 'G'];
 save(file,'G')
 file=[dir 'W' num];
 save(file,'W')
 file=[dir 'b' num];
 save(file,'b')
 file=[dir 'A' num];
 save(file,'A')
 file=[dir 'p' num];
 save(file,'p')
end
file=['x'];
filename=[dir file];
save(filename,file)
if e==1
 file='eigs_mat';
 B=[dir   file def '.fig'];
 saveas(h2,B)
 B=[dir  file def '.jpg'];
 saveas(h2,B) 
 file='eigs_mat_zoom'; 
  B=[dir   file def '.fig'];
 saveas(h21,B)
 B=[dir  file def '.jpg'];
 saveas(h21,B)
end
%break
close all
nn=5;
 iteration=500;
k1=10;
 fc=100;
 num='5';
 solve_DICCG_f(num,nn,e,iteration,dir,k1)
end

 
 