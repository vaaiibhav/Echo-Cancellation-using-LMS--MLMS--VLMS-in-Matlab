function nY = LVMakePeriodicSeq(y,N)
% LVMakePeriodicSeq([1 2 3 4],2)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
y = y(:); nY = y*([ones(1,N)]); nY = nY(:)';