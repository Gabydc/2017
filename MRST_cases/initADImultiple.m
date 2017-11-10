%function [p_ad,bhp_ad,qS_ad]=initADImultiple(p_init,qinit,n)
n=2;
foo = @(varargin)fprintf('you provided %i arguments\n',length(varargin))
%varlist2(0,0,0); 
% for i=2:n+2
%     varlist(p_ad)=varlist( 
%  
% end

s(1:4) = struct('bar',1);
foo(s.bar)


   % [varlist] = initVariablesADI(varlist2);