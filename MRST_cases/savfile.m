function savfile(dir,x,name,varargin)
nVarargs = length(varargin);
if nVarargs==0
file=[dir name];
save(file,name)
else   
file=[dir x varargin];
 save(file,x)
end