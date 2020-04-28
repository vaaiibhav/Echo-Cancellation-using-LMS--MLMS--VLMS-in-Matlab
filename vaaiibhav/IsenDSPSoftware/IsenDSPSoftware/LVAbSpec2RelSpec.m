function [Rp,As] = LVAbSpec2RelSpec(DeltaP,DeltaS)
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
%
% [Rp,As] = LVAbSpec2RelSpec(DeltaP,DeltaS)
Rp = -20*log10((1-DeltaP)/(1+DeltaP));
As = -20*log10(DeltaS/(1+DeltaP));