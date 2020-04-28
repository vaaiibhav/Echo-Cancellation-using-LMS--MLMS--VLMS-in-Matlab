function [theLocs, thePeaks] = findPeaks(theMatrix,HowMany,PeakSep)
% function [theLocs, thePeaks] = findPeaks(theMatrix,HowMany,PeakSep)
% theMatrix is a vector or matrix for which HowMany
% peak values are sought, each of which is at least PeakSep samples
% distant from the next closest peak value.  theLocs tells the locations of
% the peak values in the matrix.  Two dimensional matrices are indexed in the
% standard MATLAB manner in which the column vectors of the matrix are concatenated
% to form one long column vector, starting with the left-most column vector of the matrix.
%
% Sample Calls:
%
% [theLocs, thePeaks] = findPeaks([3 2 1; 9 7 8; 5 6 4],1,1)
%     yields theLocs=2 and thePeaks=9
% The call [theLocs, thePeaks] = findPeaks([3 2 1; 9 7 8; 5 6 4],2,2)
%     yields theLocs = [2 8] and thePeaks = [9 8]
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
aa = size(theMatrix);
theMatrix = theMatrix(1:1:aa(1)*aa(2));
nums = find(finite(theMatrix));
minVal = min(theMatrix(nums))-1;
notnums = find(~finite(theMatrix));
theMatrix(notnums) = minVal;
if HowMany>length(nums)
   HowMany=length(nums);
   Comment='Reduced number of desired peaks to the number of finite values in the test vector'
end
thePeaks(1,HowMany) = zeros;
theLocs(1,HowMany) = zeros;
for ctr = 1:1:HowMany 
[thePeaks(ctr),theLocs(ctr)] = max(theMatrix);
minLocToMinimize = max([1  (theLocs(ctr)-PeakSep)]);
maxLocToMinimize = min([length(theMatrix)  (theLocs(ctr)+ PeakSep)]);
locsToMinimize = minLocToMinimize:1:maxLocToMinimize;
theMatrix(locsToMinimize) = minVal;
end
NumNotPeaks = length(find(thePeaks==minVal));
thePeaks = thePeaks(1,1:length(thePeaks)-NumNotPeaks);
theLocs = theLocs(1,1:length(theLocs)-NumNotPeaks);
return

