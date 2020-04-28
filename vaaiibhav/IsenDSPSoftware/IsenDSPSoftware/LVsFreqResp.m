function H = LVsFreqResp(b,a,HiFreqLim,FigNo)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
    FR = (0:HiFreqLim/5000:HiFreqLim); s = j*FR;
    H = polyval(b,s)./polyval(a,s);
    figure(FigNo); subplot(311); yplot = 20*log10(abs(H)+eps);
    plot(FR,yplot); xlabel('(a) Freq, Rad/s'); ylabel('Mag, dB')
    axis([0, HiFreqLim,[min(yplot)-10], 10])
    subplot(312); plot(FR,abs(H)); xlabel('(b) Freq, Rad/s'); 
    ylabel('Mag'); axis([0 HiFreqLim 0 1.1]);
    ploty = unwrap(angle(H));
    subplot(313); plot(FR,ploty)
    xlabel('(c) Freq, Rad/s'); ylabel('Radians')
    axis([0 HiFreqLim  min(ploty)  max(ploty)])
    
 