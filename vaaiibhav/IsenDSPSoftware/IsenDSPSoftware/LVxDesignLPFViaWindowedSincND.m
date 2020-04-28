function [ImpResp,PBRipple,SBAtten, L, Status, StatusS, KaiserBeta] = LVxDesignLPFViaWindowedSincND(wp, ws, As, typeWin, LimitToTypeI)

% LVxDesignLPFViaWindowedSincND(0.2, 0.3, 60, 'blackman',0)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
if ~(strcmp(typeWin,'kaiser')==1|strcmp(typeWin,'boxcar')==1|strcmp(typeWin,'blackman')==1|...
      strcmp(typeWin,'hanning')==1|strcmp(typeWin,'hamming')==1|strcmp(typeWin,'bartlett')==1)
   fprintf('\nThe third input argument, typeWin, must be passed as any of\n')  
   fprintf('\n''kaiser'' or ''boxcar'' or ''hanning'' or ''blackman'' or\n')   
   fprintf('\n''bartlett'' or ''hamming''\n')
    Status = 0;
    StatusS = 'Invalid window type requested';
    ImpResp = 0;
    L = 0;
    PBRipple = 0;
    SBAtten = 0;
    KaiserBeta = 0;
    return
   return
end

trWidth = (ws - wp)*pi;
    
if strcmp(typeWin,'blackman')==1
     startL = ceil(12*pi/trWidth) + 1;
     minAttenDb = 74;
     KaiserBeta = [];

elseif strcmp(typeWin,'hamming')==1
     startL = ceil(8*pi/trWidth) + 1;
     minAttenDb = 53;
     KaiserBeta = [];

elseif strcmp(typeWin,'hanning')==1
     startL = ceil(8*pi/trWidth) + 1;
     minAttenDb = 44; 
     KaiserBeta = [];

elseif strcmp(typeWin,'bartlett')==1
     startL = ceil(8*pi/trWidth) + 1;
     minAttenDb = 25;
     KaiserBeta = [];

elseif strcmp(typeWin,'boxcar')==1
    minAttenDb = 21;
    startL = ceil(4*pi/trWidth) + 1;
    KaiserBeta = [];

elseif strcmp(typeWin,'kaiser')==1
    trWidth = (ws - wp)*pi;
 startL = ceil((As - 7.95)/(14.36*(trWidth/(2*pi)))) + 1;

    if As >= 50
        KaiserBeta = 0.1102*(As - 8.7);
    elseif ((As < 50) & (As >= 21))
        KaiserBeta = 0.5842*(As - 21)^0.4 + 0.07886*(As - 21);
    elseif As < 21
        KaiserBeta=0;
    end 
    minAttenDb = As;
end

if As > minAttenDb
    StatusS = 'Chosen window has its min attenuation less than As; pick another window';
    Status = 0;
    ImpResp = 0;
    L = 0;
    PBRipple = 0;
    SBAtten = 0;
    KaiserBeta = 0;
    return
end

for L = fix(startL/4)+1:1:1000
       
[ImpResp,PBRipple,SBAtten] = LVxLPFViaWindowedSincND(wp, ws, L, typeWin, KaiserBeta);

    if abs(SBAtten) >= As
            Status = 2;
            StatusS = 'L meets As requirement';
            break   
    end
    L = L + 1;
end
if abs(SBAtten) < As & L >= 1000
    Status = 1;
    StatusS = 'Reached L = 1000 without meeting As';
    ImpResp = 0;
    L = 0;
    PBRipple = 0;
    SBAtten = 0;
    KaiserBeta = 0;
    return
end
if LimitToTypeI==1
    if rem(L,2)==0
        L = L + 1;
    end
end

[ImpResp,PBRipple,SBAtten] = LVxLPFViaWindowedSincND(wp, ws, L, typeWin, KaiserBeta);
