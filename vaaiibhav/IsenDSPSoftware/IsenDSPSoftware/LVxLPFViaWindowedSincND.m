function [ImpResp,PBRipple,SBAtten] = LVxLPFViaWindowedSincND(wp, ws, L, typeWin, KaiserBeta)
% [ImpResp,PBRipple,SBAtten] = LVxLPFViaWindowedSincND(0.2, 0.3, 91,'hamming',[])
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
if ~(strcmp(typeWin,'kaiser')==1|strcmp(typeWin,'boxcar')==1|strcmp(typeWin,'blackman')==1|...
      strcmp(typeWin,'hanning')==1|strcmp(typeWin,'hamming')==1|strcmp(typeWin,'bartlett')==1)
  
   fprintf('\nThe fourth input argument, typeWin, must be passed as any of\n')  
   fprintf('\n''kaiser'' or ''boxcar'' or ''hanning'' or ''blackman'' or\n')   
   fprintf('\n''bartlett'' or ''hamming''\n')
   return
end
   
if strcmp(typeWin,'kaiser')==1
     theWin = feval('kaiser', L, KaiserBeta);
else
     theWin = feval(typeWin, L);
end

wc = (wp + ws)/2;

ImpResp = LVxTrunIdealLowpass(wc, L);

ImpResp = ImpResp.*theWin';

FFTLen = 2048;
if L >150
    FFTLen = 4096;
end

fr = fft(ImpResp,FFTLen);

fr = abs(fr);
fr = fr(1,1:fix(length(fr)/2)+1);
fr = fr/(max(fr));

PB = fr(1,1:round(wp*length(fr)));
actualPBRipple = min(PB);

PBRipple = 20*log10(actualPBRipple + eps); 

SB = fr(1,ceil(ws*length(fr)):length(fr));
actualSBAtten = max(SB);

SBAtten = 20*log10(actualSBAtten + eps);