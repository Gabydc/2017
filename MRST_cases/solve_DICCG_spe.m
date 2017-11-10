%This program gives the matrix for the solution of the Laplace problem, the
%linear system is solved with the CGF conjugate gradient function
clear all

close all
%This is the matrix that we are going to use, usually is A12
num='5';  
%The number of deflation verctors to use either 2 or 3
nn=5;
e=0;

 nxi=1;
 nyi=1;
 nx=60;
 ny=220;
  
 iteration=500;

tol=10^-5;
     close all



for k =7
    k
    close all
    dir='/home/wagm/cortes/Localdisk/Research_programs/2015_11_report/results/SPE10_60_220_85layers/';
    %dir='/dev/media/Sphinx/Doctorado_Delft/technical_report_2015/2015_11_report/results/';
 % name=['SPE10_' num2str(nx) '_' num2str(ny) ];

  %dir=[dir name '/'];

    def=['tol-' num2str(k) '_' num2str(nn) '_' num ];
text = [dir 'iter' def '.txt'];
text1 = [dir 'cond' def '.txt']; 

file=[dir 'xi' ];
load(file)
file=[dir 'A' num];
load(file)
file=[dir 'b' num];
load(file)
file=[dir 'p' num];
load(file)
file=[dir 'W' num];
load(file)
file=[dir 'G'];
load(file)
file=[dir 'l'];
load(file)
% for i=1:nn
% file=[dir 'x' num2str(i) ];
% x=load(file);
% 
% end
file=[dir 'x_' num2str(k)];
load(file)

n=length(b);

%[V,D]=defpodf_w_2(x,nn);

a=sparse(A);
l=sparse(l);
clear A
%%1
%   [V,D]=POD(x);
%  z(:,1)=x(:,5);
    s=0;
     for i=1:nn
        s=s+1;
     z(:,s)=x(:,i);
% % % %   %z(:,s)=V(:,s);
    end

%z(:,3)=x(:,3);
% %%2
% for i=1:nn
% z(:,i)=xx{i}(:,1);
% zch(:,i)=xx{i}(:,2);
% end
% clear z
%%3
% z=zeros(nx*ny,8);
% for i=1:2*nn
%     z(1+128*(i-1):n*i,i)=1;
% end
% figure
% spy(z)
% break

% clear V
% clear Vch



% 
% %  filename=[ dir '/results/qfs_1/x' num];
% % save(filename,'x')
% %  for i=1:n
% %      l(i,i)=1;
% %  end
%   for i=1:n
%      l(i,i)=a(i,i);
%   end
% l=V*inv(V'*a*V)*V';
% l=(eye(size(l))-a*l);

fileID = fopen(text1,'w');
fprintf(fileID,'\n%6s %6s %6s %6s\n\n','Method ','& ', 'Cond',' \\');
fclose(fileID);
%[x1,iter1,e1,hl1,hl11,hl111,h2,h21]=CGF(a,b,xi,iteration,tol,e,text,nn,dir,def); 

[x2,iter2,e2,hl2,h2,h21]=ICCG(a,b,xi,iteration,tol,l,e,text1,nn,dir,def,100);
         %ICCG(a,b,xi,iteration,tol,l,eig,text,nn,dir,def,fc) 
%[x3,iter3,e3,hl3,hl13,hl113,h2,h21]=DCGF(a,b,xi,iteration,tol,z,e,text1,nn,dir,def);

[x4,iter4,e4,hl4,h2,h21]=DICCG(a,b,xi,iteration,tol,z,l,e,text1,nn,dir,def,100);

 % x1(:,1)=x1;
  x2(:,1)=x2;
  %x3(:,1)=x3;
  x4(:,1)=x4;
% % filename=[dir '/results/qfs_1/xch' num];
% % save(filename,'xch')
% figure(1)
% legend([hl1,hl2,hl3,hl4],'CG','PCG','DCG','DPCG')
% ylabel('log(Error)')
% xlabel('Iteration')
% title('||x_i-x_f||_2/||x_f||_2')
figure(100)
legend([hl2,hl4],'ICCG','DICCG')
ylabel('log(||M^{-1}r^k||_2/||M^{-1}b||_2)','FontSize',15)
xlabel('Iteration','FontSize',15)
title('Convergence','FontSize',15)

% figure(6)
% legend([hl111,hl112,hl113,hl114],'CG','PCG','DCG','DPCG')
% ylabel('log(Error)')
% xlabel('Iteration')
% title('||x_i-x_f||_A/||x_f||_A')
%txt file with data
fileID = fopen(text,'w');
fprintf(fileID,'\n%6s %6s %6s %6s %6s\n\n','Method &','Iter ', '&','error', ' \\');
%fprintf(fileID,'%6s %d %6s %6.2e %6s\n','CG &',iter1,' &',max(e1), ' \\');
fprintf(fileID,'%6s %d %6s %6.2e %6s\n','PCG &',iter2,' &',max(e2), ' \\');
%fprintf(fileID,'%6s %d %6s %6.2e %6s\n','DCG &',iter3,' &',max(e3), ' \\');
fprintf(fileID,'%6s %d %6s %6.2e %6s\n','DPCG &',iter4,' &',max(e4), ' \\');
fclose(fileID);
 fprintf('\n  Method      Iteration #    error   \n');
 % fprintf('\n CG      %8d      %10.0d\n',iter1, max(e1));
  fprintf('\n ICCG %8d      %10.0d\n',iter2, max(e2));
  %fprintf('\n DCG          %8d           %1.0d\n',iter3, max(e3));
  fprintf('\n DICCG         %8d           %1.0d\n',iter4, max(e4));


% 
 figure
%h1=plotingsolution(G,W,'CG',x1,1);

h1=plotingsolution_2(G,W,'ICCG',x2,1);

%h1=plotingsolution(G,W,'DCG',x3,3);

h1=plotingsolution_2(G,W,'DICCG',x4,2);

%mkdir([dir 'solx_5_' num2str(num)])
%dir=[dir 'solx_5_' num2str(num) '/'];

mkdir([dir 'sol_1' ])
dir=[dir 'sol_1'  '/'];
% Se guarda la grafica en el directorio dir

 file='sol';
  B=[dir  file def '.fig'];
  saveas(h1,B)
% crear las carpetas para guardar los resultados



% file='eig_pod';
%  B=[dir  file def '.fig'];
%  saveas(h11,B)
%  B=[dir  file def '.jpg'];
%  saveas(h11,B)
   

% file='conv';
%  B=[dir   file def '.fig'];
%  saveas(hl4,B)
%  B=[dir   file def '.jpg'];
%  saveas(hl4,B) 
 
 file='conv_def';
 B=[dir   file def '.fig'];
 saveas(hl4,B)
 B=[dir   file def '.jpg'];
 saveas(hl4,B) 
%  file='conv_true';
%  B=[dir   file def '.fig'];
%  saveas(hl4,B)
%  B=[dir   file def '.jpg'];
%  saveas(hl4,B) 
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
 end
%   clearvars('-except','num','e','size1','iteration','tol', 'j')

