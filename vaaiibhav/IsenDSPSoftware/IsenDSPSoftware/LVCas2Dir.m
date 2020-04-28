function [b,a,k] = LVCas2Dir(Bc,Ac,Gain)
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
k = Gain; szB = size(Bc);
b = 1; for ctr = 1:1:szB(1),
b = conv(b,Bc(ctr,:)); end; 
szA = size(Ac); a = 1; 
for ctr = 1:1:szA(1),
a = conv(a,Ac(ctr,:)); end

    