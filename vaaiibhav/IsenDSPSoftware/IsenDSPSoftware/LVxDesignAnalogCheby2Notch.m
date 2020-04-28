function [Z,P,K,NetRp,NetAs] = LVxDesignAnalogCheby2Notch(OmP1,OmS1,OmS2,OmP2,Rp,As)
% Receives analog Notch band edges in radians/sec,
% plus desired maximum Rp and minimum As and returns the zeros (Z) and
% poles (P) of an analog Chebyshev Type-II Notch along with filter
% gain K and the realized values of Rp (NetRp) and As (NetAs). Also plots
% the magnitude and phase responses of both the notch filter and its
% lowpass analog prototype filter.
% Sample calls:
% [Z,P,K,NetRp,NetAs] = LVxDesignAnalogCheby2Notch(2,3,6,8,1,45)
% [Z,P,K,NetRp,NetAs] = LVxDesignAnalogCheby2Notch(1,1.2,2.8,3.9,1,60)
% [Z,P,K,NetRp,NetAs] = LVxDesignAnalogCheby2Notch(1,1.2,4,5,1,75)
% [Z,P,K,NetRp,NetAs] = LVxDesignAnalogCheby2Notch(4,5,8,10,1,40)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
    Om0Sq = OmS1*OmS2;
    OmPlp = OmP1/(Om0Sq - OmP1^2);
    OmSlp1 = OmS2/(OmS2^2 - Om0Sq);
    OmSlp2 = OmS1/(Om0Sq - OmS1^2);
    OmSlp = min([OmSlp1,OmSlp2]);
    
    [Z1,P1,K1] = LVDesignCheb2(Rp,As,OmPlp,OmSlp)
    NumPoles = length(P1)
    NumZer = length(Z1)
    Num = 1; Den = 1; 
    D = [1 0 Om0Sq]; N = [0 1 0]; 
if rem(NumPoles,2)==0  % even number of poles  
    for Ctr = 1:1:NumPoles
        Den = conv(Den,[N-P1(Ctr)*D]);
        Num = conv(Num,[N-Z1(Ctr)*D]);
    end; 
else
    for Ctr = 1:1:NumPoles
        Den = conv(Den,[N-P1(Ctr)*D]);
    end; 
    
    for Ctr = 1:1:NumZer
        Num = conv(Num,[N-Z1(Ctr)*D]);
    end;  
    Num = conv(Num,D);   
end
       
b = real(Num); a = real(Den);
s=j*[0:0.001:2*OmS2]; Hs = abs(polyval(b,s)./polyval(a,s)); 
K = 1/max(Hs); Z = roots(b); P = roots(a); 
LVsFreqRespDouble(K1*poly(Z1),poly(P1),2*OmS2,30,K*b,a);
H = LVsFreqResp(K1*b,a,1.3*OmP2,39); 
[NetAs,NetRp1,NetRp2] = LVxRealizedFiltParamNotch(H,OmP1,OmS1,OmS2,OmP2,1.3*OmP2);
NetRp = max([NetRp1,NetRp2]);