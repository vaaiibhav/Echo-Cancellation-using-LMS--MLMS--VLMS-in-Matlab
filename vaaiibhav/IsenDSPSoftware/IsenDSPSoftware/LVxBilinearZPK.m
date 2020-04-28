function [b,a] = LVxBilinearZPK(Z,P,K,Fs)
% Receives the zeros Z, poles P, and gain K of a classical IIR filter and 
% performs the Bilinear transform using Fs as the sampling rate.
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
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
    
b = real(Num); a = real(Den);

