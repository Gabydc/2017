close all
clear all
clc
pw=0;
%sz=32*2^pw;
sz=25;
%dir='/home/wagm/cortes/Localdisk/Research_programs/2016_1_report/1_22_results_temporal/';
dir='/dev/media/Sphinx/Doctorado_Delft/enero_16/results_t/';
folder=['5wells_' num2str(sz)  ];
mkdir([dir], folder)
dir = [dir folder '/'];
iter=200;
tol=10^-5;
% % Grid 
% G = cartGrid([20 20],[1 1]);
% G = computeGeometry(G);
% % load seamount
% % G = pebi(triangleGrid([x(:) y(:)], delaunay(x,y)));
% % G = computeGeometry(G);
% % 
%% Set up grid and petrophysical data
% We use a Cartesian grid of size nx-by-ny with homogeneous petrophysical
% data: permeability of 100 mD and porosity of 0.2.

nx=sz;
ny=sz;
G = cartGrid([nx, ny, 1]);
G = computeGeometry(G);

% Disable gravity
%gravity off
% Set up uniform permeability and constant porosity
%rock.perm = ones(G.cells.num, 1)*1*milli*darcy();
rock.perm=ones(G.cells.num, 1);
%inhomogeneus permeability
 lsz=round(sz*sz/8);
% for i=1:2:8
%  rock.perm(1+lsz*(i-1):lsz*i)  = repmat(10^(3)*milli*darcy(), [lsz, 1]);
% end
rock.poro = ones(G.cells.num, 1)*0.2;
% figure(10)
% plotCellData(G, rock.perm)
% break

%central well
cw=(sz*sz)-ceil((sz*sz)/2)+1;

% Grid information
N = G.faces.neighbors;
N = N(all(N ~= 0, 2), :);
nf = size(N,1);
nc = G.cells.num;

%Transmissibility
%hT = computeTrans(G, rock.perm);
%hT = computeTrans(G, struct('perm', ones(nc,1)*darcy));
hT = computeTrans(G, struct('perm', rock.perm));
cf = G.cells.faces(:,1);
T = 1 ./ accumarray(cf, 1 ./ hT, [G.faces.num, 1]);
T = T(all(N~=0,2),:);
% Operators
C = sparse([(1:nf )'; (1:nf )'], N, ones(nf,1)*[-1 1], nf, nc);
grad = @(x) C*x;
div = @(x) -C'*x;


%Constants
c = 1*10^-4;
%mu = 1*centi*poise;
mu=1;
p0 = ones(nc, 1);
%p0(r<5) = 200*atm;
% [p]=devect(G,N,nc,nf,nx,cw,T,mu);
% break
%%Snapshots
 mkdir([dir 'sol' ])
dir1=[dir 'sol'  '/'];
for i=1:5
    q1(1:4)=-0.5;
    q1(5)=1.5;
    q1(i)=0;
    if i==5
    q1(5)=2;
    end
 
    [z(:,i),A,b,h,h1]=defvect_1(p0,nc,cw,q1,div,grad,T,mu,iter,10^-7,G,nx,ny);
    file=[dir 'A' num2str(i)];
 save(file,'A')
file=[dir 'b' num2str(i)];
 save(file,'b')
 file='sol_snapshots';
B=[dir1  file num2str(i) '.fig'];
saveas(h1,B)
B=[dir1  file num2str(i) '.jpg'];
saveas(h1,B)
end


file='iter_snapshots';
B=[dir1  file  '.fig'];
saveas(h,B)
B=[dir1  file '.jpg'];
saveas(h,B)
file=[dir 'Z'];
 save(file,'z')
 close(figure(500))
%break
text=1;
nn=i;
eig=0;

q = zeros(nc, 1);
 q([1 nx cw (nc-nx+1) nc]) = [-.5 -.5  2 -.5 -.5];
presEq = @(p, p0, dt) (1/dt)*c*(p-p0) - div( (T/mu).*grad(p))+q;

 


p = initVariablesADI(p0);
[t,Tt,dt] = deal(0,10000000*day,100000*day);
s=0;
while t < Tt,
    s=s+1;
t = t + dt;
p0 = p.val;
eq = presEq(p, p0, dt);
eq([1 nx cw (nc-nx+1) nc]) = eq([1 nx cw (nc-nx+1) nc])+q([1 nx cw (nc-nx+1) nc]);
l=ichol(eq.jac{1});
[x11,h]=ICCG_0(eq.jac{1},eq.val,p0,iter,tol,l);
p.val = p.val -( x11 );
%p.val = p.val -( eq.jac{1}\eq.val );

%clf, plotCellData(G,rock.perm)
figure(1)
h2=plotCellData(G,p.val);
hold on
% plotWell(G, W);
caxis([100 200]*atm); drawnow;
 
% pause

xy(s,1)=t;
xy(s,2)=eq.val(1);
xy(s,3)=eq.val(cw);
end
file='sol_t';
B=[dir1  file  '.fig'];
saveas(h2,B)
B=[dir1  file '.jpg'];
saveas(h2,B)
figure(2); 
h3=plot(xy(:,1)/day,xy(:,2),'*',xy(:,1)/day,xy(:,3),'*') 
caxis([0 Tt]); 
ylabel('Pressure of the well','FontSize',16)
xlabel('Time (days)','FontSize',16)
legend('I','P')
file='press_wells';
B=[dir1  file  '.fig'];
saveas(h3,B)
B=[dir1  file '.jpg'];
saveas(h3,B)
hold on


% % Assemble and solve equations
% p = initVariablesADI(zeros(nc,1));
% 
% 
% %q([135 282 17]) = [-1 .5 .5];
%  % −> quarter five−spot
% eq = div(T.*grad(p))+q;
%  % equation
% eq(1) = eq(1) + p(1);
%  % make solution unique
% p = -eq.jac{1}\eq.val; % solve equation
% plotCellData(G,p);
