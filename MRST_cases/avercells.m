function [cG,crock,hg]=avercells(fine, coarse,layers,av,T)
%This program average the permeability with (1) Harmonic, (2) Arithmetic or
% (3) Harmonic-Arithmetic average.
%G = cartGrid([nx ny numel(layers)]);


    G = cartGrid([fine(1) fine(2) numel(layers)]);
    %G = cartGrid(fine, fine);
G = computeGeometry(G);
  rock = SPE10_rock(layers);
rock.perm = convertFrom(rock.perm(:,1), milli*darcy);
q = partitionUI(G,coarse);
%  cG = cartGrid(coarse,fine);
cG = cartGrid(coarse);
cG = computeGeometry(cG);
clear crock*
vol = G.cells.volumes;
%% Compute the averages
% The arithmetic and harmonic averages are straightforward. To compute the
% harmonic-arithmetic average, we  temporary partition that coincides with
% the coarse grid in along the given axial direction and with the original
% fine grid in the other directions. Then we map the averaged values back
% onto the fine grid and perform a standard arithmetic averaging.

for i=1:size(rock.perm,2)
    crock1.perm(:,i) = accumarray(q,vol.*rock.perm(:,i))./accumarray(q,vol);
    
    crock2.perm(:,i) = accumarray(q,vol)./accumarray(q,vol./rock.perm(:,i));
    
    dims = G.cartDims; dims(i)=coarse(i);
    qq = partitionUI(G, dims);
    K = accumarray(qq,vol)./accumarray(qq,vol./rock.perm(:,i));
    crock3.perm(:,i) = accumarray(q,K(qq).*vol)./accumarray(q,vol);
end
px = log10([min(rock.perm(:,1)) max(rock.perm(:,1))]);

if av==1
    crock=crock1;
else if av==2
        crock=crock2;
    else
        crock=crock3;
    end
end
figure;
hg=plotCellData(cG,log10(crock.perm(:,1)),'EdgeColor','k');
caxis(px); axis equal tight 
if T==1
title(['SPE10, grid size (' num2str(coarse(1)) ',' num2str(coarse(2)) ')'] );
else
end