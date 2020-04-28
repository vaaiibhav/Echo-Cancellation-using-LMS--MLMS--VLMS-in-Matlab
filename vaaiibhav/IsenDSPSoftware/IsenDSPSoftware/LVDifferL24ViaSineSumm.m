    function LVDifferL24ViaSineSumm
    %
    % Author: F. W. Isen
    % Copyright 2009 by Morgan & Claypool
    AkPos = -[1:1:11]*pi/12; AkLOver2 = [-12]*pi/12;
    L = 24; M = (L-1)/2; Imp = zeros(1,L); n = 0:1:L-1;
    Ak = [AkPos,AkLOver2]; limK = L/2; 
    for k = 1:1:limK
    if (k==limK); C = 1; else; C = 2; end; 
    Imp = Imp + C*Ak(k)*sin(2*pi*(n-M)*k/L); end; 
    Imp = Imp/L; figure(109); stem(Imp);
    xlabel('Sample'); ylabel('Amplitude')