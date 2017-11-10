%%If we want to save the files/graphs, it's neccesary to run this part
                
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
