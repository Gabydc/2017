function  savfig(dir,file,h,varargin)
% fprintf('Total number of inputs = %d\n',nargin);
nVarargs = length(varargin);
if nVarargs>0
 B=[dir   file varargin{1} '.fig'];
 saveas(h,B)
 B=[dir   file varargin{1} '.jpg'];
 saveas(h,B) 
else
   B=[dir   file  '.fig'];
 saveas(h,B)
 B=[dir   file '.jpg'];
 saveas(h,B) 
end
