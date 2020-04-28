function [Bbq,Abq,Gain] = LVDirToCascadeClassIIR(b,a,gain)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
    B0 = b(1); A0 = a(1); a = a/A0; b = b/B0;
    Gain = gain*B0/A0; OrdNum = length(b)-1; 
    if OrdNum == 0; Bbq = [0 0 b];
    elseif OrdNum == 1; Bbq = [0 b];
    else; [NumRts] = LVCmplxConjOrd(roots(b),10^(-10));
    LenNumRts = length(NumRts); NoPrs = fix(LenNumRts/2);
    for Ctr = 1:1:NoPrs; ind = [2*Ctr-1 2*Ctr];
    Bbq(Ctr,1:3) = real(poly(NumRts(ind))); end
    if LenNumRts-2*NoPrs > 0
    Bbq(Ctr+1,1:3) = [0 real(poly(NumRts(LenNumRts)))];
    end; end; OrdDen = length(a)-1; % denominator
    if OrdDen == 0; Abq = [0 0 a];
    elseif OrdDen == 1; Abq = [0 a];
    else; [NumRts] = LVCmplxConjOrd(roots(a),10^(-10));
    LenNumRts = length(NumRts); NoPrs = fix(LenNumRts/2);
    for Ctr = 1:1:NoPrs; ind = [2*Ctr-1 2*Ctr];
    Abq(Ctr,1:3) = real(poly(NumRts(ind))); end
    if LenNumRts-2*NoPrs > 0
    Abq(Ctr+1,1:3) = [0 real(poly(NumRts(LenNumRts)))];
    end; end; Bbq(find(abs(Bbq)<10^(-15)))=0;
    Abq(find(abs(Abq)<10^(-15)))=0; 

