function LVsFreqRespDouble(b1,a1,HiFreqLim,FigNo,b2,a2)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
    FR = (0:0.005:HiFreqLim); s = j*FR;
    H = polyval(b1,s)./polyval(a1,s);
    figure(FigNo); clf; subplot(321); 
    yplot = 20*log10(abs(H)+eps); plot(FR,yplot);
    xlabel('(a) Freq, Rad/s'); ylabel('Mag, dB');
    axis([0 inf -100 10]);
    subplot(323); plot(FR,abs(H)); 
    xlabel('(b) Freq, Rad/s'); ylabel('Mag');
    axis([0 inf 0 1.1]);
    subplot(325); plot(FR,unwrap(angle(H)))
    xlabel('(c) Freq, Rad/s'); ylabel('Radians')   
    H = polyval(b2,s)./polyval(a2,s);
    subplot(322); yplot = 20*log10(abs(H)+eps);
    plot(FR,yplot); xlabel('(d) Freq, Rad/s'); 
    ylabel('Mag, dB'); axis([0 inf -100 10]);
    subplot(324); plot(FR,abs(H)); 
    xlabel('(e) Freq, Rad/s'); ylabel('Mag');
    axis([0 inf 0 1.1]);
    subplot(326); plot(FR,unwrap(angle(H)))
    xlabel('(f) Freq, Rad/s'); ylabel('Radians')
