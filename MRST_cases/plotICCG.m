function []=plotICCG(t,ee,ne,nit)
ne=length(ee);
iter=1:ne;
if t<6
     figure(800)
     subplot(2,3,t)
     colors = [1 0 0;0 1 0; 0 1 1; 1 0 1; 0 0 1;1 0 1;0.2 0.5 0.5; 0.5 0.2 1; 1 0.2 0.5; 0.2 1 0.6;];
     h=semilogy(iter,ee,'o','Color',colors(nit,:));
     hold on

     title( ['Timestep ' num2str(t )] )
     
     %ylabel('log(||M^{-1}r^k||_2/||M^{-1}b||_2)','FontSize',16)
     %xlabel('Iteration','FontSize',16)
     else if t==25
             figure(500)
             color=[0.82 0.2 0.6];
   h=semilogy(iter,ee,'<','Color',color);
     hold on
ylabel('log(||M^{-1}r^k||_2/||M^{-1}b||_2)','FontSize',16)
xlabel('Iteration','FontSize',16)
         else if t==50
                 figure(500)
       h=semilogy(iter,ee,'pb');
     hold on
ylabel('log(||M^{-1}r^k||_2/||M^{-1}b||_2)','FontSize',16)
xlabel('Iteration','FontSize',16)
              else
         h(nit)=0;
             end 
         end    
    
     end
