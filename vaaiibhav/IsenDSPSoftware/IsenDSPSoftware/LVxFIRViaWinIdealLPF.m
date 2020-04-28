function LVxFIRViaWinIdealLPF(FiltType, BandEdgeVec, typeWin, As, Rp)
% function LVxFIRViaWinIdealLPF(FiltType, BandEdgeVec, typeWin, As, Rp)
% FiltType: 1 = Lowpass, 2 = Highpass, 3 = Bandpass, 4 = Notch
%
% BandEdgeVec must have 2 or 4 values; it must have length 2 when defining ws and wp for an LPF or
% an HPF; it must have length 4 for a BPF or Notch filter
%
% typeWin is one of the following functions in string format, i.e., surrounded by single quotes
% rectwin, kaiser, blackman, hanning, hamming, bartlett
%
% As is minimum desired stopband attenuation in dB, as a positive number
% such as 51, 67, etc.
%
% Typical calls:
%
% Lowpass filters
%  LVxFIRViaWinIdealLPF(1, [0.2 0.3], 'kaiser', 40, 0.1)
%  LVxFIRViaWinIdealLPF(1, [0.2 0.3], 'kaiser', 60, 0.1)
%  LVxFIRViaWinIdealLPF(1, [0.2 0.3], 'kaiser', 80, 0.1)
%  LVxFIRViaWinIdealLPF(1, [0.2 0.3], 'boxcar', 21, 0.1)
%  LVxFIRViaWinIdealLPF(1, [0.2 0.3], 'hanning', 44, 0.1)
%  LVxFIRViaWinIdealLPF(1, [0.2 0.3], 'hamming', 53, 0.1)
%  LVxFIRViaWinIdealLPF(1, [0.2 0.3], 'blackman', 74, 0.1)
% 
% Highpass filters
%  LVxFIRViaWinIdealLPF(2, [0.5 0.6], 'kaiser', 40, 0.1)
%  LVxFIRViaWinIdealLPF(2, [0.5 0.6], 'kaiser', 60, 0.1)
%  LVxFIRViaWinIdealLPF(2, [0.5 0.6], 'kaiser', 80, 0.1)
%  LVxFIRViaWinIdealLPF(2, [0.5 0.6], 'boxcar', 21, 0.1)
%  LVxFIRViaWinIdealLPF(2, [0.5 0.6], 'hanning', 44, 0.1)
%  LVxFIRViaWinIdealLPF(2, [0.5 0.6], 'hamming', 53, 0.1)
%  LVxFIRViaWinIdealLPF(2, [0.5 0.6], 'blackman', 74, 0.1)
%
% Bandpass filters
%  LVxFIRViaWinIdealLPF(3, [0.2 0.3 0.5 0.6], 'kaiser', 40, 0.1)
%  LVxFIRViaWinIdealLPF(3, [0.2 0.3 0.5 0.6], 'kaiser', 60, 0.1)
%  LVxFIRViaWinIdealLPF(3, [0.2 0.3 0.5 0.6], 'kaiser', 80, 0.1)
%  LVxFIRViaWinIdealLPF(3, [0.2 0.3 0.5 0.6], 'boxcar', 21, 0.1)
%  LVxFIRViaWinIdealLPF(3, [0.2 0.3 0.5 0.6], 'hanning', 44, 0.1)
%  LVxFIRViaWinIdealLPF(3, [0.2 0.3 0.5 0.6], 'hamming', 53, 0.1)
%  LVxFIRViaWinIdealLPF(3, [0.2 0.3 0.5 0.6], 'blackman', 74, 0.1)
%
% Bandstop filters
%  LVxFIRViaWinIdealLPF(4, [0.2 0.3 0.5 0.6], 'kaiser', 40, 0.1)
%  LVxFIRViaWinIdealLPF(4, [0.2 0.3 0.5 0.6], 'kaiser', 60, 0.1)
%  LVxFIRViaWinIdealLPF(4, [0.2 0.3 0.5 0.6], 'kaiser', 80, 0.1)
%  LVxFIRViaWinIdealLPF(4, [0.2 0.3 0.5 0.6], 'boxcar', 21, 0.1)
%  LVxFIRViaWinIdealLPF(4, [0.2 0.3 0.5 0.6], 'hanning', 44, 0.1)
%  LVxFIRViaWinIdealLPF(4, [0.2 0.3 0.5 0.6], 'hamming', 53, 0.1)
%  LVxFIRViaWinIdealLPF(4, [0.2 0.3 0.5 0.6], 'blackman', 74, 0.1)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

if ~(strcmp(typeWin,'kaiser')==1|strcmp(typeWin,'boxcar')==1|strcmp(typeWin,'blackman')==1|...
   strcmp(typeWin,'hanning')==1|strcmp(typeWin,'hamming')==1|strcmp(typeWin,'bartlett')==1)
   fprintf('\nThe third input argument, typeWin, must be passed as any of\n')  
   fprintf('\n''kaiser'' or ''boxcar'' or ''hanning'' or ''blackman'' or\n')   
   fprintf('\n''bartlett'' or ''hamming''\n')
   return
end

BandEdgeVec = sort(BandEdgeVec);
lenBEV = length(BandEdgeVec);

if lenBEV==2
    if ~(FiltType==1|FiltType==2)
        error('Length of BandEdgeVec was 2, but FiltType was not 1 or 2 (LPF or HPF)')
    end
elseif lenBEV==4
    if ~(FiltType==3|FiltType==4)
        error('Length of BandEdgeVec was 4, but FiltType was not 3 or 4 (BPF or Bandstop)')
    end 
end

if (FiltType==1|FiltType==2)
    wp = BandEdgeVec(1);
    ws = BandEdgeVec(2);
        if FiltType==1
        [ImpResp,PBRipple,SBAtten, L, Status, StatusS, KaiserBeta] = LVxDesignLPFViaWindowedSincND(wp, ws, As, typeWin, 0);
            if Status==0
            error('Initial L inadequate to meet As, pick different window (such as Kaiser), or reduce As')
            end
        else
        [ImpResp,PBRipple,SBAtten, L, Status, StatusS, KaiserBeta] = LVxDesignLPFViaWindowedSincND(wp, ws, As, typeWin, 1);
            if Status==0
            error('Initial L inadequate to meet As, pick different window (such as Kaiser), or reduce As')
            end
        end
     if FiltType==2 % convert to Highpass
        UnityL = zeros(1,L);
        M = (L-1)/2  + 1;
        UnityL(1,M) = 1;
        ImpResp = UnityL - ImpResp;
        dummy = wp;
        wp = ws;
        ws = dummy;
      end
end

% A Bandstop filter will be synthesized as the sum of the impulse responses of
% a lowpass filter using wp and ws as the lowest two values of BandEdgeVec,
% and a highpass filter using the highest two values of BandEdgeVec

% A Bandpass filter will be synthesized as an LPF impulse response using
% the two largest values of BandEdgeVec, minus an LPF impulse response
% using the two lowest values of BandEdgeVec. 

% For either filter type, design both LPF's first and pick the larger
% resulting value of L, then redesign the shorter LPF if necessary using
% the function
% [ImpResp,PBRipple,SBAtten] = LPFViaWindowedSincND(wp, ws, L, typeWin, KaiserBeta)

if FiltType==3 % BPF 
    wp3 = BandEdgeVec(3);
    ws4 = BandEdgeVec(4);
   [ImpResp,PBRipple,SBAtten, L1, Status, StatusS, KaiserBeta] = LVxDesignLPFViaWindowedSincND(wp3, ws4, As, typeWin, 0);
        if Status==0
            error('Initial L inadequate to meet As, pick different window (such as Kaiser), or reduce As')
        end
    wp1 = BandEdgeVec(1);
    ws2 = BandEdgeVec(2);
   [ImpResp,PBRipple,SBAtten, L2, Status, StatusS, KaiserBeta] = LVxDesignLPFViaWindowedSincND(wp1, ws2, As, typeWin, 0);
        if Status==0
            error('Initial L inadequate to meet As, pick different window (such as Kaiser), or reduce As')
        end
   L = max([L1  L2]);
   actualSBAtten = 0;
while abs(actualSBAtten) < As
   [WideImp,PBRipple,SBAtten] = LVxLPFViaWindowedSincND(wp3, ws4, L, typeWin, KaiserBeta);
   [NarrowImp,PBRipple,SBAtten] = LVxLPFViaWindowedSincND(wp1, ws2, L, typeWin, KaiserBeta);   
   ImpResp = WideImp - NarrowImp;  % the net Bandpass Impulse Response
      [actualPBRipple, actualSBAtten] = BandpassEval(ImpResp,BandEdgeVec);
     if abs(actualSBAtten) < As
        L = L + 1;
     end
end  
end

if FiltType == 4 % Bandstop % limit to odd length
    wp1 = BandEdgeVec(1);
    ws2 = BandEdgeVec(2); % get L = Lnarr for the narrow LPF
   [ImpResp,PBRipple,SBAtten, Lnarr, Status, StatusS, KaiserBeta] = LVxDesignLPFViaWindowedSincND(wp1, ws2, As, typeWin, 1);
        if Status==0
            error('Initial L inadequate to meet As, pick different window (such as Kaiser), or reduce As')
        end
    wp3 = BandEdgeVec(3);
    ws4 = BandEdgeVec(4); % get L = Lwide for the wide LPF
   [ImpResp,PBRipple,SBAtten, Lwide, Status, StatusS, KaiserBeta] = LVxDesignLPFViaWindowedSincND(wp3, ws4, As, typeWin, 1);
        if Status==0
            error('Initial L inadequate to meet As, pick different window (such as Kaiser), or reduce As')
        end
L = max([Lnarr  Lwide]);
actualSBAtten = 0;
while abs(actualSBAtten) < As
    [WideImp,PBRipple,SBAtten] = LVxLPFViaWindowedSincND(wp3, ws4, L, typeWin, KaiserBeta);
    % convert to HPF impulse response
            UnityL = zeros(1,L);
            M = (L-1)/2  + 1;
            UnityL(1,M) = 1;
            HPFImpResp = UnityL - WideImp;
    % get the narrow LPF
    [NarrowImp,PBRipple,SBAtten] = LVxLPFViaWindowedSincND(wp1, ws2, L, typeWin, KaiserBeta);  
    ImpResp = HPFImpResp + NarrowImp; % net Bandstop filter impulse response
    
     [actualPBRipple, actualSBAtten] = BandstopEval(ImpResp,BandEdgeVec);
     if abs(actualSBAtten) < As
         if rem(L,2)==0 % L even
             L = L + 1;
         else
         L = L + 2;
         end
     end
end
end

figure(10951)
clf

% Evaluate passband ripple and stopband attenuation and display
FFTLen = 2048;
if L >150
    FFTLen = 4096;
end
fr = fft(ImpResp,FFTLen);

fr = abs(fr);
fr = fr(1,1:fix(length(fr)/2)+1);
fr = fr/(max(fr));

if FiltType==1
PB = fr(1,1:round(wp*length(fr)));
PBRipple = min(PB);
actualPBRipple = 20*log10(PBRipple + eps); %------------FiltType 1--PB

SB = fr(1,round(ws*length(fr)):length(fr));
SBAtten = max(SB);
actualSBAtten = 20*log10(SBAtten + eps);   %------------FiltType 1--SB
end

if FiltType==2
wp = BandEdgeVec(2);
PB = fr(1,round(wp*length(fr)):length(fr));
PBRipple = min(PB);
actualPBRipple = 20*log10(PBRipple + eps); %------------FiltType 2--PB

ws = BandEdgeVec(1);
SB = fr(1,1:round(ws*length(fr)));
SBAtten = max(SB);
actualSBAtten = 20*log10(SBAtten + eps);   %------------FiltType 2--SB
end

if FiltType==3 % one passband, two stopbands
    [actualPBRipple, actualSBAtten] = BandpassEval(ImpResp,BandEdgeVec);
end


if FiltType==4 % two passbands, one stopband
[actualPBRipple, actualSBAtten] = BandstopEval(ImpResp,BandEdgeVec);
end

actualPBRipple =  abs(actualPBRipple);
actualSBAtten = abs(actualSBAtten);

plot( [0:1:length(fr)-1]/length(fr), 20*log10(fr + eps) )

if strcmp(typeWin,'kaiser')==1
    xlabel(['Normalized Frequency (L = ',num2str(L),');',' ( Window = kaiser ({\beta} = ',num2str(KaiserBeta),') )'])
else
xlabel(['Normalized Frequency (L = ',num2str(L),');',' ( Window = ',typeWin,' )'])
end
ylabel('20log10(Mag(DTFT))')

if FiltType==1
text(0.02,35,['Max Passband Ripple = ',num2str(actualPBRipple,2),' dB'])
text(0.52,35,['Min Stopband Atten = ',num2str(actualSBAtten,2),' dB'])
text(0.02,25,['Design Rp = ',num2str(Rp,2),' dB'])
text(0.52,25,['Design As = ',num2str(As,2),' dB'])
text(0.15,10,['wp = ',num2str(wp),'\pi',' Radians'])
text(0.65,10,['ws = ',num2str(ws),'\pi',' Radians'])
elseif FiltType==2
text(0.52,35,['Max Passband Ripple = ',num2str(actualPBRipple,2),' dB'])
text(0.02,35,['Min Stopband Atten = ',num2str(actualSBAtten,2),' dB'])
text(0.52,25,['Design Rp = ',num2str(Rp,2),' dB'])
text(0.02,25,['Design As = ',num2str(As,2),' dB'])
text(0.65,10,['wp = ',num2str(wp),'\pi',' Radians'])
text(0.15,10,['ws = ',num2str(ws),'\pi',' Radians'])    
else
text(0.02,35,['Max Passband Ripple = ',num2str(actualPBRipple,2),' dB'])
text(0.52,35,['Min Stopband Atten = ',num2str(actualSBAtten,2),' dB'])
text(0.02,25,['Design Rp = ',num2str(Rp,2),' dB'])
text(0.52,25,['Design As = ',num2str(As,2),' dB'])
BE = BandEdgeVec;
strBE = [num2str(BE(1)),'  ',num2str(BE(2)),'  ',num2str(BE(3)),'  ',num2str(BE(4))];
text(0.12,10,['Band Edges = [ ',strBE,' ] times \pi',' Radians'])
end
axis([0 1 -max([120 abs(As)]) 45])

function [actualPBRipple, actualSBAtten] = BandstopEval(ImpResp,BandEdgeVec)
L = length(ImpResp);
FFTLen = 2048;
if L >150
    FFTLen = 4096;
end
fr = fft(ImpResp,FFTLen);

fr = abs(fr);
fr = fr(1,1:fix(length(fr)/2)+1);
fr = fr/(max(fr));

%if FiltType==4 % two passbands, one stopband
% first passband
wplo = 1;
wphi = BandEdgeVec(1);
PBlo = fr(1,1:round(wphi*length(fr)) );
actualPBRipple1 = min(PBlo);

% second passband
wphi = BandEdgeVec(4);
PBhi = fr(1,round(wphi*length(fr)):length(fr) );
actualPBRipple2 = min(PBhi);

PBRipple = min([actualPBRipple1 actualPBRipple2 ]);
actualPBRipple = 20*log10(PBRipple + eps);%-----------------FiltType 4--PB

% stopband
wslo = BandEdgeVec(2);
wshi = BandEdgeVec(3);

SB = fr(1,round(wslo*length(fr)):round(wshi*length(fr)) );
SBRipple = max(SB);

actualSBAtten = round(20*log10(SBRipple + eps)); %----------------FiltType 4--SB

function [actualPBRipple, actualSBAtten] = BandpassEval(ImpResp,BandEdgeVec)
L = length(ImpResp);
FFTLen = 2048;
if L >150
    FFTLen = 4096;
end
fr = fft(ImpResp,FFTLen);

fr = abs(fr);
fr = fr(1,1:fix(length(fr)/2)+1);
fr = fr/(max(fr));
% pb
wplo = BandEdgeVec(2);
wphi = BandEdgeVec(3);
PB = fr(1,round(wplo*length(fr)):round(wphi*length(fr)) );
PBRipple = min(PB);
actualPBRipple = 20*log10(PBRipple + eps); %--------------FiltType 3--PB

% first stopband
wshi = BandEdgeVec(1);
SBlo = fr(1,1:round(wshi*length(fr)) );
SBAtten1 = max(SBlo);

% second stopband
wslo = BandEdgeVec(4); 
SBhi = fr(1,round(wslo*length(fr)):length(fr) );
SBAtten2 = max(SBhi);

netSB = max([SBAtten1 SBAtten2]);
actualSBAtten = 20*log10(netSB + eps); %-------------------FiltType 3--SB
