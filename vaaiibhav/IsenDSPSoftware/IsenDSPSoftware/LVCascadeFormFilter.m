function [y] = LVCascadeFormFilter(Bc,Ac,Gain,x)
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
szA = size(Ac);
for ctr = 1:1:szA(1)
a = Ac(ctr,:); b = Bc(ctr,:);
x = filter(b,a,x);
end
y = Gain*x;