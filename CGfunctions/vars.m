%% Set up grid and petrophysical data
% We use a Cartesian grid of size nx-by-ny with homogeneous petrophysical
% data: permeability of 100 mD and porosity of 0.2.

%Select a linear solver 1.AGMG, 2.GMRES, 3.PCG_IC, 4.DPCG_IC, other Backslash
lsolver=3;
% Select transport solver 1. Explicit, 2. Implicit
tsolver=2;
% We define the size of the reservoir and the inicial points of the grid
% sz is the number of cells, is the same for the two dimensions
sz=64;

 
% If we use wells, we can define here the position of the wells
% w1=10*2^pw;
% w2=2*w1;
nxi=1;
nyi=1;
nx=sz;
ny=sz;
nz=1;
layers = 1:85;
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

% Number of deflation vectors
dv=10;
podv{1}=[1:10];
podv{2}=[6:10];