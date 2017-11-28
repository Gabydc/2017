function TOF_F(state, G, rock,varargin)
opt = struct( 'reverse', false, ...
    'times' , false, ...
    'dir', [], ...
    'plots', false, ...
    'nf', '0', ...
    'wells', [],...
    'bc', [],...
    'src', []);
opt = merge_options(opt, varargin{:});
nf1 = opt.nf;
dir = opt.dir;
W  = opt.wells;
src = opt.src;
bc = opt.bc;

if(~isempty(W))
    T  = computeTimeOfFlight(state, G, rock, 'W',W);
     name = '-W';
else if(~isempty(src))
        T  = computeTimeOfFlight(state, G, rock, 'src',src);
        name = '-src';
    else if(~isempty(bc))
            T  = computeTimeOfFlight(state, G, rock, 'bc',bc);
             name = '-bc';
        end
    end
end

nf = str2num(nf1)-1;
nf1 = nf;
cmin = min(T);
cmax = max(T);

if(opt.plots)
    nf = nf + 1;
    f(nf)=figure(nf);
    file{nf} = ['Time-of-flight' name];
    plotCellData(G, T, 'edgecolor','k','edgealpha',0.05);
    title([file{nf} ' is: ' num2str(cmax)]);
    caxis([cmin,cmax]);axis equal tight;colormap jet
    colorbar       
    if(opt.times)
        disp([file{nf} ' is: ' num2str(cmax)])
     end
    if(opt.reverse)
        if(~isempty(W))
            RT  = computeTimeOfFlight(state, G, rock, 'W',W,'reverse',true);
           
        else if(~isempty(src))
                RT  = computeTimeOfFlight(state, G, rock, 'src',src,'reverse',true);
           
            else if(~isempty(bc))
                RT  = computeTimeOfFlight(state, G, rock, 'bc',bc,'reverse',true);
               
                end
            end
        end
        
        
        cmin = min(RT);
        cmax = max(RT);
      
        nf = nf + 1;
        f(nf)=figure(nf);
        file{nf} = ['Reverse-Time-of-flight' name];
        plotCellData(G, RT, 'edgecolor','k','edgealpha',0.05);
        title([file{nf} ' is: ' num2str(cmax)]);
        caxis([cmin,cmax]);axis equal tight;colormap jet
        colorbar
        if(opt.times)
        disp([file{nf} ' is: ' num2str(cmax)])
        end
        
        cmin = min(T+RT);
        cmax = max(T+RT);
        nf = nf + 1;
        f(nf)=figure(nf);
        file{nf} = ['Total-time' name];
        plotCellData(G, T+RT, 'edgecolor','k','edgealpha',0.05);
        title([file{nf} ' is: ' num2str(cmax)]);
        caxis([cmin,cmax]);axis equal tight;colormap jet
        colorbar
          if(opt.times)
        disp([file{nf} ' is: ' num2str(cmax)])
        end
    end
  
    for i = nf1+1: nf
        f(i) = figure(i);
        savefigures(f(i), file{i}, dir)
    end
    
    
    
end


