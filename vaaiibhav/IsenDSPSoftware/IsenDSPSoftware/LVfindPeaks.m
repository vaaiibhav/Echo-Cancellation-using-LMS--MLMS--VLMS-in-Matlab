function [Locs, Pks] = LVfindPeaks(Mat,QuanPks,PeakSep)
% function [Locs, Pks] = LVfindPeaks(Mat,QuanPks,PeakSep)
% Mat is a vector or matrix for which QuanPks peak values are sought, each 
% of which is at least PeakSep samples distant from the next closest peak 
% value. Locs tells the locations of the peak values in the matrix. To
% make adjacent samples of Mat eligible to be peaks, use PeakSep = 0.
% Sample call:
% [Locs, Pks] = LVfindPeaks([0:1:7],8,0)
% [Locs, Pks] = LVfindPeaks([0:1:7],8,1)
% [Locs, Pks] = LVfindPeaks((cos(2*pi*3*[0:1:63]/64)),4,8)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

aa = size(Mat); 
Mat = Mat(:);  
nums = find(isfinite(Mat));
minVal = min(Mat)-1;
notnums = find(isnan(Mat));
Mat(notnums) = minVal;
    if QuanPks > length(nums), 
        QuanPks = length(nums);
     Comment = 'Num peaks reduced'; 
    end
Pks = zeros(1,QuanPks);
Locs = zeros(1,QuanPks);

for ctr = 1:1:QuanPks 
  [aa,bb] =  max(Mat);
  Pks(1,ctr) = aa;
  Locs(1,ctr) = bb;
  minLocToMin = max([1, (Locs(ctr)-PeakSep)]);
  maxLocToMin = min([length(Mat),(Locs(ctr)+ PeakSep)]);
  locsToMin = minLocToMin:1:maxLocToMin;
  Mat(locsToMin) = minVal; 
end;

NumNotPeaks = length(find(Pks==minVal));
Pks = Pks(1,1:length(Pks)-NumNotPeaks);
Locs = Locs(1,1:length(Locs)-NumNotPeaks);