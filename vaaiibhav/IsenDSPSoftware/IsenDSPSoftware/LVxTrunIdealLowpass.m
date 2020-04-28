function b = LVxTrunIdealLowpass(fwc, L)
% function b = LVxTrunIdealLowpass(fwc, L)
%
% Returns a truncated ideal lowpass impulse response
%
% fwc is a fraction of pi radians; for cutoff of radian frequency pi/2, use fwc = 0.5;
% L is length of desired Impulse Response
%
% Typical call:
%
% b = LVxTrunIdealLowpass(0.5, 51)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
M = (L-1)/2; wc = fwc*pi; n = 0:1:L-1;
b = sin(wc*(n - M + eps))./(pi*(n - M + eps)); 
