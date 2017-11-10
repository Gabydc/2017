clear all
close all
%condition to compute the eigenvalues, ife=1, the eigenvalues will be
%computed
e=0;
%The directory where we want to store the matrices, for each size a new
%directory will be created
dir='/home/wagm/cortes/Localdisk/Results/16_03/POD/';
%  dir='/dev/media/Sphinx/Doctorado_Delft/Research_programs/2016_1_report/';
  mrstModule add spe10 coarsegrid;
 nxi=1;
 nyi=1;
 nx=60;
 ny=220;
 nxf=60;
 nyf=220;
 nz=1;
 nzf=1;
 layers = 2;
 fine = [nxf nyf nzf];
 coarse = [nx ny nz];
  name=['SPE10_' num2str(nx) '_' num2str(ny) ];
  mkdir([dir], name )
  dir=[dir name '/'];
iteration=200;
  [G,rock,hg]=avercells(fine, coarse,layers,3,0);
  savfig(dir,'Grid',hg)
for k=1:2:7
    pause
    tol=10^-k;
%num is the number of snapshot
 for num=1:5
     %It is neccesary to create a random vector as initial guess and save
     %it for the iteration
    if num==1
        xi=rand(nx*ny,1);
%         savfile(dir,xi,'xi')
         file=[dir 'xi'];
         save(file,'xi')
    end 
    
    
%     figure
%     plotCellData(G,log(permupscaled));
    
%% Define well locations
wtype    = {'bhp', 'bhp', 'bhp', 'bhp', 'bhp'};
wtarget  = [-1,   -1,   -1,  -1,  3] .* barsa();
wrad     = [0.125, 0.125, 0.125, 0.125, 0.125] .* meter;
wloc     = [  nxi,   nx,     nxi,   nx,   nx/2;
              nyi,   ny,     ny,   nyi,   ny/2];
wname    = {'P1', 'P2', 'P3', 'P4', 'I1'};
sgn      = [ -1 ,  -1 ,  -1 ,  -1 ,   1 ];


wtarget(num)=0;
if num==5
   wtarget(num)=4*barsa();
end



W = [];


for w = 1 : numel(wtype),
   W = verticalWell(W, G, rock, wloc(1,w), wloc(2,w), [], ...
                    'Type', wtype{w}, 'Val', wtarget(w), ...
                    'Radius', wrad(w), 'Name', wname{w}, ...
                    'Sign', sgn(w), 'InnerProduct', 'ip_tpf');
end

% figure
% h = plotCellData(G,log(rock.perm)); 
% h = plotWell(G, W);
% axis equal tight; colormap(jet(128));
% view(60,30)
% file='perm_layer_';
% B=[dir   file  '.fig'];
% saveas(h(1),B)
% B=[dir   file  '.jpg'];
% saveas(h(1),B)
 





%Boundary conditions
%  bc  = pside([], G, 'YMin',10.*barsa());
%  bc  = pside(bc, G, 'YMax',0.*barsa());
sol = initState(G, W, 0);
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
% is present)
gravity reset off
fluid = initSingleFluid('mu' , 1*centi*poise, ...
                        'rho', 1014*kilogram/meter^3);
%display(fluid)

psolve = @(state) incompTPFA(state,G, hT, fluid, 'wells', W,'MatrixOutput',true);

 sol= psolve(sol);
 nu=num2str(num);
 
 
 
 

%h=plotWell(G, W);
view(0,90)
p=sol.pressure;
 A=sol.A(1:G.cells.num,1:G.cells.num);
 b=sol.rhs(1:G.cells.num);
  def=['A' nu];
text = [dir 'cond.txt'];
text1 = [dir 'eigs' def  '.txt'];


if num==1
l=ichol(A);
 file=[dir 'l'];
 save(file,'l')
else
    file=[dir 'l'];
load(file)
end


[x11,hc0,hc01,hc02]=ICCG_0(A,b,xi,iteration,tol,l,e); 


if e==1
     file='eigs_a';
  B=[dir  file '.fig'];
  saveas(hc01(1),B)
   B=[dir  file '.jpg'];
  saveas(hc01(1),B)
end


 x(:,num)=x11;
 x1(:,num)=sol.pressure;
 res(:,num)=x(:,num)-x1(:,num);
dir1=[dir 'sol_snap/'];
mkdir([dir1])

figure(num+200)
[ht]=plotingsolution(G,W,'direct', sol.pressure,1) ;
colorbar
title('backslash')
savfig(dir1,'sol_bs',ht,nu)
figure(num+300)
 [ha]=plotingsolution(G,W,'ICCG',x(:,num),1);
 colorbar
title('ICCG')
 [ha]=plotingsolution(G,W,'bs-ICGCG',res(:,num),2);
 colorbar
title('bs-ICCG')
savfig(dir1,'sol_app',ha,nu)


 file=[dir 'G'];
 save(file,'G')
 file=[dir 'W' nu];
 save(file,'W')
 file=[dir 'b' nu];
 save(file,'b')
 file=[dir 'A' nu];
 save(file,'A')
 file=[dir 'p' nu];
 save(file,'p')
 
if e==1
 [Va,Da,h2,h21]=eigvs(A,text,text1,num,6);  
 savfig(dir,'eigs',h2,def)
 savfig(dir,'eigs_zoom',h21,def)
 savfig(dir,'Da',h2,def)
%  file='eigs';
%  B=[dir   file def '.fig'];
%  saveas(h2,B)
%  B=[dir   file def '.jpg'];
%  saveas(h2,B) 
%   file='eigs_zoom';
%  B=[dir   file def '.fig'];
%  saveas(h21,B)
%  B=[dir   file def '.jpg'];
%  saveas(h21,B) 
 file=[dir 'D' num];
 save(file,'Da')
end

 end
 %%
 file=['x'];
filename=[dir file];
save(filename,file)
file=['x1'];
filename=[dir file];
save(filename,file) 

 close all
nn=4;
 iteration=500;
k1=7;
 fc=100;
 defv='4';
 solve_DICCG_f_pod(nu,nn,e,iteration,dir,k1);
 end
