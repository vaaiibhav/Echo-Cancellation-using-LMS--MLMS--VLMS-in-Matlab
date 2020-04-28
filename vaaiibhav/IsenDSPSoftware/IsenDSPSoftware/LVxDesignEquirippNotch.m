function [actRp,actAs,WF] = LVxDesignEquirippNotch(Rp,As,wp1,ws1,ws2,wp2)
% function [actRp,actAs,WF] = LVxDesignEquirippNotch(Rp,As,wp1,ws1,ws2,wp2)
% Rp is the maximum design passband ripple in positive dB and As is the minimum design stopband
% attenuation in dB, wp1, ws1,ws2, and wp2 specify the 1st passband, 1st stopband, 2nd stopband, & 2nd
% passband frequencies, respectively. actRp and actAs are the realized values of Rp and As, 
% respectively, and bis the vector of designed filter coefficients.
% [actRp,actAs,WF] = LVDesignEquirippNotch(0.2,60,0.4,0.45,0.65,0.7) 
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
    BndLm = [0,wp1,ws1,ws2,wp2,1];
    Rfac= 10^(-Rp/20); DeltaP = (1-Rfac)/(1+Rfac);
    DeltaS = (1+DeltaP)*10^(-As/20);
    deltaF = abs(wp1 - ws1)/2;
    compL = (-20*log10(sqrt(DeltaP*DeltaS)) - 13)/(14.6*deltaF) + 1;
    L = round(compL), Ord = L - 2; LenFFT = 2^15;
    SBAtten = 0; PBR = 100;
    if rem(Ord,2)==0 % order must be even, i.e. length must be odd
    else
        Ord = Ord + 1;
    end
    while (SBAtten < As)|(PBR > Rp); Ord = Ord + 2;
    WF = firpm(Ord,BndLm,[1,1,0,0,1,1],[DeltaS/DeltaP,1,DeltaS/DeltaP]);
    fr = abs(fft(WF,LenFFT)); 
    fr = fr(1,1:LenFFT/2+1)/max(fr);
    Lfr = length(fr); 
    
    PB1 = fr(1,1:floor(wp1*Lfr));  
    PB2 = fr(1,ceil(wp2*Lfr):Lfr);
    Rp1 = -20*log10(min(PB1)+eps); 
    Rp2 = -20*log10(min(PB2)+eps);
    PBR = max(Rp1,Rp2);
    
    SB = fr(1,ceil(ws1*Lfr):floor(ws2*Lfr)); 
    SBAtten = -20*log10(max(SB)+eps);
    end; 
    
    Fin_L = Ord + 1, actRp = PBR; actAs = SBAtten; figure(8); clf;
    plot([0:1:LenFFT/2]/(LenFFT/2), 20*log10(fr+eps)); 
    xlabel('Normalized Frequency (Units of \pi)')
    ylabel('Magnitude, dB'); axis([0,1,-(As+20),5])