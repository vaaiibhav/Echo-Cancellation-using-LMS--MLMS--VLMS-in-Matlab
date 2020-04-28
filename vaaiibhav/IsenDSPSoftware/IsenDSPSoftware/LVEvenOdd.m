function [xe,xo] = LVEvenOdd(x)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
xe = (x + fliplr(x))/2; xo = (x - fliplr(x))/2;