%% Time-of-flight
% To study how heterogeneity affects  flow patterns and define natural
% time-lines in the reservoir, it is common to study the so-called
% time-of-flight (TOF), i.e., the time it takes an imaginary particle
% released at an inflow boundary or at a perforation of an injector to
% reach a given point in the reservoir. Time-of-flight are usually
% associated with streamline methods, but can also be computed from linear
% steady-state transport equations on the form,
%
% $$v\cdot \nabla \tau = \phi. $$
%
% In this example, we will show how to compute time-of-flight from a given
% flux field.
clear all
close all
clc
mrstModule add incomp diagnostics streamlines

%% Setup model
% As our model, we will use a logically Cartesian grid in which the
% node perturbed so that grid lines form curves and not lines. For
% simplicity, we assume unit permeability and porosity and use a set of
% source and sink terms to emulate a quater five-point setup.
[nx,ny,nz] = deal(35,35,1);
per = 0;
cp = 0;
%Create the grid
% G = twister(G);
% G = computeGeometry(G);
G = cartGrid([nx, ny, nz], [nx, ny, nz]);
G = computeGeometry(G);



%% Forward time-of-flight
clf

%rock = makeRock(G, 1, 1);


% Set up uniform permeability and constant porosity
rock.perm = ones(G.cells.num, 1)*10*milli*darcy();
rock.poro = ones(G.cells.num, 1)*0.2;

% Create layers of permeability
nlay = 5;
v = [];
for i = 6:2*nlay:ny
    for j = 0:nlay-1
        if i+j < ny+1
            v = [v j+i];
        end
    end
end


% [I] = Sub2ind([1:nx],1:ny,v,nx,ny,nz);
[I] = Sub2ind([1:nx],v,1:nz,nx,ny,nz);
rock.perm(I) = 10*10^(-per)*milli*darcy();
figure(199)
file{1} = ['Permeability'];
plotCellData(G, rock.perm/(milli*darcy())); colorbar
%return
axis equal tight off
title('Permeability field [mD]')
% return
if nz > 1
    view(10,20)
else
    view(0,90)
end




pc_form = 'nonwetting';
[muw, muo] = deal( 1, 10);
[rhow, rhoo] = deal(1000, 700);
[krw, kro] = deal(2, 2);
cap_scale = 10;
props = constantProperties([   muw,  muo] .* centi*poise, ...
    [rhow, rhoo] .* kilogram/meter^3);
%[kr, pc]  = tabulatedSatFunc([x, x.^krw, y.^kro, y.*cap_scale*barsa]);
if cp==0
    fluid = initSimpleFluid('mu' , [   muw,    muo] .* centi*poise     , ...
        'rho', [rhow, rhoo] .* kilogram/meter^3, ...
        'n'  , [   krw,    kro]);
else
    
    fluid = struct('properties', props                  , ...
        'saturation', @(x, varargin)    x.s  , ...
        'relperm'   , kr                     , ...
        'pc'        , @(x, varargin) pc(x.s));
    xDummy   = initState(G, [], [0, 1]);
    xDummy.s = linspace(0, 1, numel(xDummy.s))'; ...
        pc = convertTo(fluid.pc(xDummy), barsa);
    nf = nf + 1;
    hcp=figure(nf);
    plot(xDummy.s, pc);
    xlabel('s_w'); ylabel('pc [bar]');
    title('Capillary pressure curve')
    file{nf} = ['Capillary pressure'];
end
pv = poreVolume(G, rock);
%injRate = -sum(pv)/(50*day);
injRate = -0.4*meter^3/day;
bc = fluxside([], G, 'xmin', -injRate, 'sat', [1, 0]);
bc = pside(bc, G, 'xmax', 0*barsa, 'sat', [0, 1]);
src = addSource([], 1, sum(poreVolume(G,rock)), 'sat', 1);
src = addSource(src, G.cells.num, -sum(poreVolume(G, rock)), 'sat', 1);





%% Solve pressure equation
% Compute transmissibilities and solve pressure equation
trans  = computeTrans(G, rock);
xr = incompTPFA(initResSol(G, 100), G, trans, fluid, 'src', src);
xb = incompTPFA(initResSol(G, 100), G, trans, fluid, 'bc', bc);
dir ='/mnt/sda2/cortes/Programs/2017/MRST_programs/Injection_B/Results/';
TOF_F(xr, G, rock,'reverse', true, ...
    'times' , true, ...
    'plots', true, ...
    'nf', '1', ...
    'src', src, 'dir', dir);

TOF_F(xb, G, rock,'reverse', true, ...
    'times' , true, ...
    'plots', true, ...
    'nf', '10', ...
    'bc', bc, 'dir', dir);




%% Copyright notice

% <html>
% <p><font size="-1">
% Copyright 2009-2017 SINTEF ICT, Applied Mathematics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
