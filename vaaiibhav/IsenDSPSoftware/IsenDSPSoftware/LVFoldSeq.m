function [xFold,nFold] = LVFoldSeq(x,n)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
xFold = fliplr(x), nFold = -fliplr(n)