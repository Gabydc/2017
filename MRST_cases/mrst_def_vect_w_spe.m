%This program gives the matrix for the solution of the Laplace problem, the
%linear system is solved with the CGF conjugate gradient function
clear all
close all
dir='/home/wagm/cortes/Localdisk/Research_programs/2015_11_report/results_2/';
 %dir='/dev/media/Sphinx/Doctorado_Delft/technical_report_2015/2015_11_report/results/';
 nxi=1;
 nyi=1;
 nx=30;
 ny=110;
  name=['SPE10_' num2str(nx) '_' num2str(ny) ];

  dir=[dir name '/'];
iteration=500;

nn=5;
for t=1:2:9
tol=10^-t;
fileA=[dir 'A5'];
load(fileA); 

for i=1:nn
num=num2str(i); 
fileb=[dir 'b' num];
load(fileb)
filep=[dir 'p' num]; 
load(filep)
files=[dir 'W' num];
load(files)
fileG=[dir 'G'];
load(fileG)

% [V,D] = eigs(A);
% figure(2)
% h1=subplot(2,2,1);
% hold on
% plot(diag(D),'*r'); 
% title('Eigenvalues');
% ylabel('Value ')
% xlabel('Eigenvalue')
if i==1
xi(1:length(b),1)=rand;
file=['xi'];
filename=[dir file];
save(filename,file)
l=ichol(A);
file=['l'];
filename=[dir file];
save(filename,file)
else 
   file=[dir 'xi'];
load(file)
file=[dir 'l'];
load(file)
end
[x11]=ICCG_0(A,b,xi,iteration,tol,l);
 x(:,i)=x11;
 figure(i+10*t)
 [h1]=plotingsolution(G,W,'CG',x(:,i),2);
end
save([dir 'x_' num2str(t) ],'x')

end