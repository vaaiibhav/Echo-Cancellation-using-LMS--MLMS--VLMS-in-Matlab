function LVCheby1toCheby2(Ep,N)
% function LVCheby1toCheby2(Ep,N)
%
% Converts a Chebyshev Type-I Magnitude-Squared function to a Chebyshev
% Type-II Magnitude-Squared function
% Ep is the parameter epsilon, and N is filter order
%
% Sample Call:
% LVCheby1toCheby2(0.5,5)
% 
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
inc = 0.005; 
xLo = 0:inc:1;
xHi = 1+inc:inc:3; 
TnLo = cos(N*acos(xLo));
TnHi = cosh(N*acosh(xHi)); 
T = [TnLo TnHi];
MagHSq = 1./(1 + Ep^2*(T.^2));
figure(333)

subplot(211); 
    xplot = [xLo, xHi];
    plot(xplot,MagHSq); 
    xlabel('Freq, Rad/s')
    ylabel('Mag Squared')
    % Convert to Type-II
    xLo = xLo(find(~(xLo==0))); 
    TnLo = cos(N*acos(1./xLo));
    xHi = xHi(find(~(xHi==0))); 
    TnHi = cosh(N*acosh(1./xHi));
    T = [TnLo TnHi]; 
    MagHSq = 1 - 1./(1 + Ep^2*(T.^2));
subplot(212); 
    xplot = [xLo, xHi];
    plot(xplot,MagHSq); 
    xlabel('Freq, Rad/s')
    ylabel('Mag Squared')