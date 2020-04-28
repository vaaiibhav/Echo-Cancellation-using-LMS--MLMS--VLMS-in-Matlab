function [DeltaP,DeltaS] = LVRelSpec2AbSpec(Rp,As)
% [DeltaP,DeltaS] = LVRelSpec2AbSpec(0.5,60)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
Rfac= 10^(-Rp/20); DeltaP = (1-Rfac)/(1+Rfac);
DeltaS = (1+DeltaP)*10^(-As/20);