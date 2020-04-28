function LVRealFiltfromCCPoles(PoleMag,PoleAng)
% LVRealFiltfromCCPoles(0.95,pi/4)
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
    Pole = PoleMag*exp(j*PoleAng); cPole = conj(Pole);
    rcoeffs = conv([1, -(Pole)],[1, -(cPole)]); SR = 1024; 
    t = 0:1/(SR-1):1; x = chirp(t,0,1,SR/2); 
    y = filter([1],[rcoeffs],x);
    figure(8); plot(y); 
    xlabel('Sample'); ylabel('Amplitude')
    axis([0, length(y), 1.1*min(y), 1.1*max(y)])