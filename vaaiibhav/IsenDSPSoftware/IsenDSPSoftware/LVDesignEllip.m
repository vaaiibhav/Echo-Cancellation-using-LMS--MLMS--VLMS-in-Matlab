function [Z,P,K] = LVDesignEllip(Rp,As,OmgP,OmgS) 
% [Z,P,K] = LVDesignEllip(1.25,50,0.5,0.6)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
E = (10^(Rp/10)-1)^0.5;
A=10^(As/20); OmgC = OmgP;
k = OmgP/OmgS; k1 = E/sqrt(A^2-1);
[K1k, K2k] = ellipke([k, (1-k^2)].^2);
[K1k1, K2k1] = ellipke([k1, (1-k1^2)].^2);
N = ceil(K1k(1)*K1k1(2)/(K1k1(1)*K1k(2)));   
[Z,P,K] = LVellip(N,Rp,As,OmgC);
    
   

    
