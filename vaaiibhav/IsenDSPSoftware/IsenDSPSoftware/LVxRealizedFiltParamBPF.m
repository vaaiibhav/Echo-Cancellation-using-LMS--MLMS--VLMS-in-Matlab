function [NetRp,NetAs1,NetAs2] = LVxRealizedFiltParamBPF(H,OmS1,OmP1,OmP2,OmS2,HiFreqLim)
% Receives the four band edge design frequencies for a bandpass filter, OmS1,OmP1,PmP2,and OmS2,
% the computed complex frequency response H from frequency 0 up to
% HiFreqLim, and returns the realized value of Rp (NetRp) for the passband and the
% realized values of As (NetAs1 and NetAs2) for each stopband.
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
    Lfr = length(H);
    absH = abs(H);
    H = H/max(absH);
    
    mFr = HiFreqLim;
    
    sb1End = floor(OmS1/mFr*Lfr); 
    sb1 = absH(1:sb1End); 
    
    pbStart = ceil((OmP1/mFr)*Lfr); 
    pbEnd = floor((OmP2/mFr)*Lfr);
    pb = absH(pbStart:pbEnd); 
    
    sb2Start = ceil((OmS2/mFr)*Lfr);  
    sb2 = absH(sb2Start:Lfr);
   
    mnpb = min(pb);
    NetRp = 20*log10(1/mnpb);
    
    maxsb1 = max(sb1); 
    maxsb2 = max(sb2); 

    NetAs1 = 20*log10(1/maxsb1);
    NetAs2 = 20*log10(1/maxsb2);