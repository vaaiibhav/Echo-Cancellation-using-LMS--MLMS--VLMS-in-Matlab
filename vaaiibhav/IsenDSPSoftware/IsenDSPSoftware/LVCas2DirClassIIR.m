function [b,a,k] = LVCas2DirClassIIR(Bbq,Abq,Gain)
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
k = Gain; szB = size(Bbq);
b = 1; for ctr = 1:1:szB(1),
b = conv(b,Bbq(ctr,:)); end; 
szA = size(Abq); a = 1; 
for ctr = 1:1:szA(1),
a = conv(a,Abq(ctr,:)); end