function [Va,Da,h2,h21]=eigvs(a,text,text1, nam,ev)
n=size(a,2);

[Va,Da] = eigs(a,n);
 Da=diag(abs(Da));
 conda=condest(a,2)
 condadeff=max(Da)/min(Da)
 
fileID = fopen(text,'a');
fprintf(fileID,'%6s %6s %6.2e %6s\n',nam,' &',condadeff, ' \\');
fclose(fileID);
fileID = fopen(text1,'w');
fprintf(fileID,'%6.2e \n', Da);
fclose(fileID);
figure(3)

Da=log(abs(Da));
h2=plot(Da,'*');
axis('tight')
t=['Eigenvalues ' nam];
title(t)
ylabel('log scale')
xlabel('Eigenvalue')
hold on
figure(4)
title(t)
subplot(2,1,1)

h21=plot(Da,'*');
if ev>1
axis([1, ev, min(Da(1:ev)),max(Da(1:ev)) ])
end
ylabel('log scale')
xlabel('Eigenvalue ')
title('First Eigenvalues (Zoom)')
hold on
subplot(2,1,2)

h21=plot(Da,'*');
axis([n-ev, n, min(Da(n-ev:n)),max(Da(n-ev:n)) ])
ylabel('log scale')
xlabel('Eigenvalue')
title('Last Eigenvalues (Zoom)')
hold on



end