function [NetRp,NetAs] = LVzRealizedFiltParamLPF(Hz,OmgP,OmgS,HiFreqLim)
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
    Lfr = length(Hz); mFr = HiFreqLim;
    Lenpb = round(OmgP/mFr*Lfr); pb = Hz(1:Lenpb); 
    mnpb = min(abs(pb)); NetRp = 20*log10(1/mnpb);
    sbStrt = (OmgS/mFr); sb = Hz(round(sbStrt*Lfr):Lfr);
    maxsb = max(abs(sb)); NetAs = 20*log10(1/maxsb);
    
    