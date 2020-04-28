function [y] = LVxFilterParallelForm(Bp,Ap,Cp,x)
% function [y] = LVxFilterParallelForm(Bp,Ap,Cp,x)
% Receives a set of Parallel Form coefficients and a sequence x
% and filters x in a Parallel Form structure, delivering the filtered
% output
% as y
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
szA = size(Ap);
y = 0;
for ctr = 1:1:szA(1)
    y = y + filter(Bp(ctr,:),Ap(ctr,:),x);
end
if ~(isempty(Cp))
    y = y + filter(Cp,1,x);
end


    