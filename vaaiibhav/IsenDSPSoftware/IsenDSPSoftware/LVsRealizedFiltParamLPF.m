function [NetRp,NetAs] = LVsRealizedFiltParamLPF(H,OmgP,OmgS,HiFreqLim)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
    Lfr = length(H); mFr = HiFreqLim;
    Lenpb = round(OmgP/mFr*Lfr); pb = H(1:Lenpb); 
    mnpb = min(abs(pb)); NetRp = 20*log10(1/mnpb);
    sbStrt = (OmgS/mFr); sb = H(round(sbStrt*Lfr):Lfr);
    maxsb = max(abs(sb)); NetAs = 20*log10(1/maxsb);