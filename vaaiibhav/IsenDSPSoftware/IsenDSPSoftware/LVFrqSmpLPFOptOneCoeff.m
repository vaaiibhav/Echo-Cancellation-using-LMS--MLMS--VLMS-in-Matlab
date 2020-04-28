function  [NormSampFreq,limK,BestAk] = LVFrqSmpLPFOptOneCoeff(wp,ws,L,Thigh,Tlow,Dec,pausedur)
% [NormSampFreq,limK,Ak] = LVFrqSmpLPFOptOneCoeff(0.5,0.6,40,0.45,0.35,-0.01,[1])
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
ctr = 0; 
    ftsz1 = 16;
    noComp = ceil((Thigh-Tlow)/abs(Dec)) 
    TSB = zeros(noComp,3);
    LenFFT = 2048; 
    if rem(L,2)==0; 
        limK = L/2-1;
    else;
    limK = (L-1)/2;
    end; 
    n = 0:1:L-1;
    M=(L-1)/2;
    NormSampFreq = [0:1:limK]/(L/2)

    wc = (wp + ws)/2;
    passInd = find(NormSampFreq <= wp);
    maxPBInd = max(passInd);
  
    stopInd = find(NormSampFreq >= ws);
    szPBInd = size(passInd);
    szSBInd = size(stopInd);
    lenPB = max([szPBInd(1),szPBInd(2)]);
    lenSB = max([szSBInd(1),szSBInd(2)]);

    for T = Thigh:Dec:Tlow; 

    Ak = [ones(1,lenPB), T, zeros(1,lenSB)];
    LA = length(Ak);
    WF = (cos(((n-M)')*[0:1:LA-1]*2*pi/L))*([Ak(1) 2*Ak(2:LA) ]' );
    WF = WF'; 
    Imp = WF/L;
    figure(3)
  
    subplot(211); 
    stem([0:1:length(Imp)-1],Imp); 
    xlabel('(a) n'); 
    ylabel('Amplitude'); 
    
    fr = abs(fft(Imp,LenFFT)); 
    fr = fr(1,1:LenFFT/2+1); 
    fr = fr/(max(abs(fr))); 
    Lfr = length(fr); 
    SB = fr(1,round(ws*(Lfr-1)):Lfr); 
%    PB = fr(1,1:round((wp/pi)*(Lfr-1)));
    PB = fr(1,1:round(wp*(Lfr-1)));
    SBAt = -20*log10(max(SB) +eps);
    PBR = -20*log10(min(PB)); 
    ctr = ctr + 1; 
    TSB(ctr,1) =T; 
    TSB(ctr,2)=SBAt; 
    TSB(ctr,3)=PBR;
    
 subplot(212)
    plot([0:1:Lfr-1]/(Lfr-1),20*log10(fr+eps)); 
    strSBA = num2str(SBAt);
    
    %hold on
    %line([0,1],[-SBAt,-SBAt]);
    strStats = ['T = ',num2str(T),' and As = ',strSBA];
    text(0.05,-SBAt+8,strStats)
    xlabel(['(b) Norm Freq, T = ',num2str(T),' and As = ',strSBA])
    ylabel('Magnitude'); 
    axis([0 1 -70 10]);
    if isempty(pausedur)
        pause(1)
    else
    pause(pausedur);
    end
    end; 
    sc = TSB(:,2); 
    bestSB = max(sc),
    bestSBind = find(sc==bestSB); 
    bestT = TSB(bestSBind,1),
    finalPBRipple = TSB(bestSBind,3)
    BestAk = [ones(1,lenPB), bestT, zeros(1,lenSB)];