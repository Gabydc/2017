% We set up the reservoir properties, the domain and the boundaries/wells

SPE = true;
Heter = false;
Plots = false;
nf = 0;
Wells = true;


steps = 30;

%% Set up model
% Square domain with homogeneous rock properties, no flow across the
% boundaries, no gravity, injection in the southeast and production in the
% northwest corners, both operating at fixed bottom-hole pressure
[nxi, nyi, nzi] = deal(1,1,1);

domain = [250,250,200];
per =1;
if Heter == true
    [nx, ny, nz] = deal(35,35,5);
    cartDim = [nx ny nz];
    G      = computeGeometry(cartGrid(cartDim,domain));
    rock   = makeRock(G, 100*milli*darcy, 0.2);
    
    %We construct a vector that contains the number of the cells that have
    %different permeability
    cperl = 5;  % Number of cells in one layer
    ilc = 6;    % Initial layer cell
    v = [];
    for i = ilc:2*cperl:ny
        for j = 0:cperl-1
            if i+j < ny+1
                v = [v j+i];
            end
        end
    end
    
    % [I] = Sub2ind([1:nx],1:ny,v,nx,ny,nz);
    [I] = Sub2ind([1:nx],v,1:nz,nx,ny,nz);
    rock.perm(I) = 10*10^(-per)*milli*darcy();
    
end
if(SPE)
    layers =1;
    
    [nx, ny, nz] = deal(60,220,1);
    cartDim = [nx, ny, numel(layers)];
    domain = cartDim .* [20, 20, 2]*ft;
    rock = SPE10_rock(layers);
    is_pos             = rock.poro > 0;
    rock.poro(~is_pos) = min(rock.poro(is_pos));
    G      = computeGeometry(cartGrid(cartDim,domain));
    % max(rock.perm)/min(rock.perm)
    %The permeability is
    %perm = rock.perm(:,1);
    %perm = reshape(perm,[60 220 1]);
    %permupscaled = sampleFromBox(G,perm);
    %rock.perm=permupscaled;
    %max(rock.perm)
    %min(rock.perm)
    %contrast=max(rock.perm)/min(rock.perm)
end
pv     = poreVolume(G, rock);
if(Plots)
    nf = nf + 1;
    f(nf) = figure(nf);
    figure(nf)
    file{nf} = 'Permeability';
    h=plotCellData(G, log(rock.perm(:,1)),'LineStyle','none');
  %  h=plotCellData(G, log(rock.perm(:,1)));
    axis equal tight off
end

if(Wells)
    %% Wells
    nw = 5;
    units_w = barsa;
    well = zeros(nw,1);
    well(1:nw) = [1100 275 275 275 275];
    wtype    = {'bhp', 'bhp', 'bhp', 'bhp', 'bhp'};
    wtarget  = [well(1),   well(2),   well(3),   well(4), well(5)].*units_w;
    wrad     = [0.125, 0.125, 0.125, 0.125, 0.125] .* meter;
    wloc     =  {Sub2ind(nx,ny,1:nz,nx,ny,nz), ...
        Sub2ind(nxi,nyi,1:nz,nx,ny,nz),Sub2ind(nxi,ny,1:nz,nx,ny,nz), ...
        Sub2ind(nx,nyi,1:nz,nx,ny,nz),Sub2ind(ceil(nx/2),ceil(ny/2),1:nz,nx,ny,nz)};
    wname    = {'P1', 'P2', 'P3', 'P4','I1'};
    Sign      = [ -1 ,  -1 ,  -1 ,  -1 , 1 ];
    Comp_iI= [1 0];
    Comp_iP= [0 1];
    
    % Extra well
%     wtype{6} = 'bhp';
%     wtarget(6) = well(6)*units_w;
%     wrad(6) = 0.125*meter;
%     wloc{6} = Sub2ind(ceil((nx+2)/2),ceil((ny+2)/2),1:nz,nx,ny,nz);
%     wname{6} = 'I2';
%     Sign(6) = 1;
    
    for i = 1 : numel(wtype)
        if Sign(i) == -1
            wcomp{i} = Comp_iP;
        else
            wcomp{i} = Comp_iI;
        end
    end
   % nw=2;
    W = [];
    for w = 1 : nw
        W = addWell(W, G, rock, wloc{w}, ...
            'Type', wtype{w}, 'Val', wtarget(w), ...
            'Radius', wrad(w), 'Name', wname{w}, ...
            'Comp_i', wcomp{w}, ...
            'Sign', Sign(w));
    end
    
    
    
    if(Plots)
        file{nf} = 'Permeability_wells';
        for i = 1 : nw
            plotWell(G, W(i),'color', 'r');
        end
        plotWell(G, W(1), 'color', 'b');
%        plotWell(G, W(6), 'color', 'b');
        
        if nz > 1
            view(10,20)
        else
            view(0,90)
        end
        hold off
    end
    
  
    %Properties to plot wells
    pmark =['o' '+' '*' 'x' 'd' 's' ];
    pcol ={'r' [0 0.7 0] 'b'  [0 0.7 0.7] [0.7 0 0.7] [0.5 0.5 0.7] };
    
else if(Boundary)
        %injRate = -sum(pv)/(500*day);
        %injRate = -0.4*meter^3/day;
        bc = fluxside([], G, 'xmin', -injRate, 'sat', [1, 0]);
        bc = pside(bc, G, 'xmax', 0*barsa, 'sat', [0, 1]);
    end
end
