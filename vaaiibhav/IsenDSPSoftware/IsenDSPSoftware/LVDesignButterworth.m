function [Z,P,K] = LVDesignButterworth(OmgP,OmgS,Rp,As)
% [Z,P,K] = LVDesignButterworth(0.4,0.7,0.2,40)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
    num = 10^(Rp/10)-1; denom = 10^(As/10)-1;
    N = ceil(log10(num/denom)/(2*log10(OmgP/OmgS)));
    OmcP = OmgP/((10^(Rp/10)-1)^(0.5/N)) 
    OmcS = OmgS/((10^(As/10)-1)^(0.5/N)) 
    OmgC = (OmcP + OmcS)/2; k = 0:1:2*N-1; 
    P = OmgC*exp(j*(pi/(2*N))*(2*k+N+1));
    P = P(find(real(P)<0)); Z = [];
    K = OmgC^N;
