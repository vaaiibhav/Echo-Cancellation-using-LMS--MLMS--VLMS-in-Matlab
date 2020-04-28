function [NetAs,NetRp1,NetRp2] = LVxRealizedFiltParamNotch(H,OmP1,OmS1,OmS2,OmP2,HiFreqLim)
% Receives the four band edge design frequencies for a notch filter, OmP1,OmS1,PmS2,and OmP2,
% the computed complex frequency response H from frequency 0 up to
% HiFreqLim, and returns the realized values of Rp (NetRp1 and NetRp2) for each passband and the
% realized value of As (NetAs) for the stopband.
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
    Lfr = length(H); absH = abs(H);
    H = H/max(absH); mFr = HiFreqLim;  
    pb1End = floor(OmP1/mFr*Lfr); 
    pb1 = absH(1:pb1End); 
    sbStart = ceil((OmS1/mFr)*Lfr); 
    sbEnd = floor((OmS2/mFr)*Lfr);
    sb = absH(sbStart:sbEnd);   
    pb2Start = ceil((OmP2/mFr)*Lfr);  
    pb2 = absH(pb2Start:Lfr);  
    mnsb = max(sb); NetAs = 20*log10(1/mnsb);    
    minpb1 = min(pb1); minpb2 = min(pb2); 
    NetRp1 = 20*log10(1/minpb1);
    NetRp2 = 20*log10(1/minpb2);