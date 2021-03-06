function [b,a,G,NetRp,NetAs] = LVxDesignDigCheby2Notch(Rp,As,wp1,ws1,ws2,wp2)
% Receives digital Notch filter band edges in normalized frequency (units of pi),
% plus desired maximum Rp and minimum As and returns the b (numerator) and
% a (denominator) coefficients of a digital Cheby2 Notch filter along with filter
% gain G and the realized values of Rp (NetRp) and As (NetAs).
% Sample calls:
% [b,a,G,NetRp,NetAs] = LVxDesignDigCheby2Notch(1,45,0.4,0.475,0.65,0.775)
% [b,a,G,NetRp,NetAs] = LVxDesignDigCheby2Notch(1,60,0.45,0.55,0.7,0.83)
% [b,a,G,NetRp,NetAs] = LVxDesignDigCheby2Notch(1,75,0.4,0.475,0.65,0.775)
% [b,a,G,NetRp,NetAs] = LVxDesignDigCheby2Notch(1,75,0.1,0.12,0.4,0.5)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

Fs = 1;
T = 1/Fs; 
OmP1 = (2/T)*tan(wp1*pi/2);
OmS1 = (2/T)*tan(ws1*pi/2);
OmS2 = (2/T)*tan(ws2*pi/2);
OmP2 = (2/T)*tan(wp2*pi/2);

[Z,P,K,NetRp,NetAs] = LVxDesignAnalogCheby2Notch(OmP1,OmS1,OmS2,OmP2,Rp,As);

NumPoles = length(P);
NumZeros = length(Z);

Num = K; Den = 1; 

if NumPoles==NumZeros
 
Den = 1; 
for Ctr = 1:1:length(P); 
    N = [(2*Fs-Z(Ctr)), -(2*Fs+Z(Ctr))];
    Num = conv(Num,N); 
    D = [(2*Fs-P(Ctr)), -(2*Fs+P(Ctr))];
    Den = conv(Den,D); 
end; 
    
elseif (NumPoles - NumZeros) == 1
    for Ctr = 1:1:NumPoles;  
        D = [(2*Fs-P(Ctr)), -(2*Fs+P(Ctr))];
        Den = conv(Den,D); 
    end;  
    for Ctr = 1:1:NumZeros;  
        N = [(2*Fs - Z(Ctr)), -(2*Fs + Z(Ctr))];
        Num = conv(Num,N); 
    end;  
    Num = conv(Num,[1,1]);
else
    Comment = 'Number of poles and zeros abnormal, cannot compute transfer function'
    return
end
    
zB = real(Num); zA = real(Den);
b0 = 1/zB(1); a0 = 1/zA(1);
a = zA*a0; b = zB*b0;
zzs = 0:pi/2000:pi;
fr = polyval(b,exp(j*zzs))./polyval(a,exp(j*zzs));
G = 1/max(abs(fr));
Hz = LVzFr(G*b,a,2000,48);
[NetAs,NetRp1,NetRp2] = LVxRealizedFiltParamNotch(Hz,wp1,ws1,ws2,wp2,1);
NetRp = max([NetRp1,NetRp2]);


