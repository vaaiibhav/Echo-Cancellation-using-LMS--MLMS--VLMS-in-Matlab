function WF = LVxFilterViaCosineFormula(Type,Bin0,PosBins)
% WF = LVxFilterViaCosineFormula(Type,Bin0,PosBins)
% WF = LVxFilterViaCosineFormula(2,[1],[1,1,1,1,0,0,0,0,0])
% A Type I filter is obtained by passing Type as 1
% A Type II filter is obtained by passing Type as 2
% Ak is broken into Bin0 = Ak[0] and all other positive bins, passed as
% PosBins
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
    Ak = [Bin0,PosBins]; if Type==1; L = 2*(length(Ak)-1)+1; 
    else; L = 2*length(Ak); end; M = (L-1)/2; 
    WF = zeros(1,L); n = 0:1:L-1;
    if rem(L,2)==0; LimK = L/2 - 1; else;
    LimK = (L-1)/2; end; LA = length(Ak);
    WF = ((2*(cos((n-M)'*(2*pi*(1:1:LimK)/L)))*Ak(1,2:LA)' + Ak(1))/L)';
    fr = fft(WF,2048); fr=fr(1,1:1025); 
    figure(59); subplot(221); stem(0:1:L-1,WF); xlabel('(a) Sample'); 
    xvec = [0:1:1024]/1024; ylabel(['Amplitude']); 
    subplot(222); plot(xvec,abs(fr)); grid on;
    xlabel('(b) Frequency, Units of \pi'); ylabel(['Magnitude']);
    subplot(223); plot(xvec,unwrap(angle(fr)));
    xlabel('(c) Frequency, Units of \pi'); ylabel(['Radians']);
    subplot(224); plot(xvec,20*log10(abs(fr)+eps)); grid on;
    xlabel('(d) Frequency, Units of \pi'); ylabel(['Magnitude, dB']); 
    axis([0 1 -100 20]);