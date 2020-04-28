% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
% this script is not meant to be called; it contains clean copies of m-code
% from the text. To run any of the example m-code segments below, copy the
% code, paste it into the Command line, and press Return.
return
% Page 101-----------------------------------------------------------------
b = [1,0,7.6088]; a = [1,0.9783,1.2435,0.5265];
[Bbq,Abq,Gain]=LVDirToCascadeClassIIR(b,a,1)
% Page 132-----------------------------------------------------------------
[b,a] = butter(5,pi,'s')
H= LVsFreqResp(b,a,2*pi,7);
% Page 133-----------------------------------------------------------------
LVButterMagSqCurves
% Page 135-----------------------------------------------------------------
LVButterPoles(3,2)
% Page 136-----------------------------------------------------------------
LVButterFR(3,1)
%--------------------------------------------------------------------------
LVButterFRViaPoly(3,1,4)
% Page 138-----------------------------------------------------------------
[z,p,k] = buttap(3)
%--------------------------------------------------------------------------
b = 27; a = poly([(-1.5 + j*2.598),(-1.5 - j*2.598),-3]);
H= LVsFreqResp(b,a,20,13);
%--------------------------------------------------------------------------
LVButterPolesAndSysFcn(3,3)
% Page 139-----------------------------------------------------------------
[Z,P,K] = LVDesignButterworth(0.4,0.7,0.2,40)
% Page 140-----------------------------------------------------------------
Rp = 0.2; As = 40; OmgP = 0.4; OmgS = 0.7;
[Z,P,K] = LVDesignButterworth(OmgP,OmgS,Rp,As)
H= LVsFreqResp(K*poly(Z),poly(P),4*OmgS,14)
[NetRp,NetAs] = LVsRealizedFiltParamLPF(H,OmgP,...
OmgS,4*OmgS)
[Bbq,Abq,Gain] = LVDirToCascadeClassIIR(poly(Z),poly(P),K)
% Page 141-----------------------------------------------------------------
LVCheby1MagSquared(0.4,5)
% Page 143-----------------------------------------------------------------
LVCheby1(5,1,0.5)
% Page 145-----------------------------------------------------------------
p=LVCheby1Poles2ndMethod(5,2,0.5)
%--------------------------------------------------------------------------
[z,p,k] = cheb1ap(2,0.2)
%--------------------------------------------------------------------------
LVCheby1PolesAndSysFcn(2,4,0.2)
% Page 146-----------------------------------------------------------------
[Z,P,K] = LVDesignCheby1Filter(0.5,40,0.5,0.65)
% Page 147-----------------------------------------------------------------
Rp = 0.5; As = 40; OmgP = 0.5; OmgS = 0.65;
[Z,P,K] = LVDesignCheby1Filter(Rp,As,OmgP,OmgS)
H= LVsFreqResp(K*poly(Z),poly(P),2*OmgS,15);
[NetRp,NetAs] = LVsRealizedFiltParamLPF(H,OmgP,...
OmgS,2*OmgS)
[Bbq,Abq,Gain] = LVDirToCascadeClassIIR(poly(Z),poly(P),K)
% Page 148-----------------------------------------------------------------
LVCheby1toCheby2(0.5,5)
% Page 149-----------------------------------------------------------------
[z,p,k] = cheb2ap(3,40)
% Page 150-----------------------------------------------------------------
LVCheb2(3,40,2)
% Page 151-----------------------------------------------------------------
[Z,P,K] = LVDesignCheb2(0.2,40,0.9,1)
% Page 152-----------------------------------------------------------------
Rp = 0.2; As = 40; OmgP = 0.9; OmgS = 1;
[Z,P,K] = LVDesignCheb2(Rp,As,OmgP,OmgS)
H= LVsFreqResp(K*poly(Z),poly(P),2*OmgS,17);
[NetRp,NetAs] = LVsRealizedFiltParamLPF(H,OmgP,OmgS,2*OmgS)
[Bbq,Abq,Gain] = LVDirToCascadeClassIIR(poly(Z),poly(P),K)
% Page 153-----------------------------------------------------------------
LVellip(5,0.2,40,2)
% Page 154-----------------------------------------------------------------
[k,e] = ellipke(0.5^2)
% Page 155-----------------------------------------------------------------
[Z,P,K] = LVDesignEllip(1.25,50,0.5,0.6)
%--------------------------------------------------------------------------
Rp = 1.25; As = 50; OmgP = 0.5; OmgS = 0.6;
[Z,P,K] = LVDesignEllip(Rp,As,OmgP,OmgS)
H= LVsFreqResp(K*poly(Z),poly(P),3*OmgS,19);
[NetRp,NetAs] = LVsRealizedFiltParamLPF(H,OmgP,...
OmgS,3*OmgS)
[Bbq,Abq,Gain] = LVDirToCascadeClassIIR(poly(Z),poly(P),K)
% Page 157-----------------------------------------------------------------
[z,p,k] = buttap(2); p = 2*p; b1=4;
a1 = poly(p); p = [2.121*(-1 + j), 2.121*(-1 - j)];
b2=9; a2 = poly(p);
LVsFreqRespDouble(b1,a1,6,22,b2,a2)
% Page 159-----------------------------------------------------------------
[z,p,k] = buttap(2); p = 2*p; b1=4;
a1 = poly(p); P = roots([1,4.2426,9]);
z = [0 0]; b2 = poly(z); a2 = poly(P);
LVsFreqRespDouble(b1,a1,6,25,b2,a2)
% Page 162-----------------------------------------------------------------
[Z,P,K] = LVCheb1Lpf2Hpf(7,1,3,3)
% Page 164-----------------------------------------------------------------
[Z,P,K] = LVDesignCheby1BPF(4,5,8,10,1,40)
% Pages 164-165------------------------------------------------------------
OmS1=4; OmP1=5; OmP2=8; OmS2=10; Rp=1; As=40;
[Z,P,K] = LVDesignCheby1BPF(OmS1,OmP1,OmP2,OmS2,Rp,As)
[Bbq,Abq,Gain] = LVDirToCascadeClassIIR(poly(Z),poly(P),K)
H= LVsFreqResp(K*poly(Z),poly(P),3*OmS2,19);
[NetRp,NetAs1,NetAs2] = LVxRealizedFiltParamBPF(H,...
OmS1,OmP1,OmP2,OmS2,3*OmS2)
% Page 166-----------------------------------------------------------------
[Z,P,K] = LVDesignCheby1Notch(4,5,8,10,1,40)

[Z,P,K] = LVDesignCheby1Notch(4,5,8,10,1,40)
% Page 170-----------------------------------------------------------------
LVDesignDigFiltViaImpInv(0.4*pi,0.5*pi,1,40,1,1)
% Page 171-----------------------------------------------------------------
LVDesignDigFiltViaImpInv(0.5*pi,0.7*pi,0.5,40,1,1)
% -------------------------------------------------------------------------
LVDesignDigFiltViaImpInv(0.5*pi,0.65*pi,0.5,40,1,1)
% -------------------------------------------------------------------------
LVDesignDigFiltViaImpInv(0.5*pi,0.62*pi,0.5,40,1,1)
% -------------------------------------------------------------------------
LVDesignDigFiltViaImpInv(0.5*pi,0.7*pi,0.5,40,2,1)
% -------------------------------------------------------------------------
LVDesignDigFiltViaImpInv(0.5*pi,0.525*pi,0.5,40,2,1)
% Page 173-----------------------------------------------------------------
LVDesignDigFiltViaImpInv(0.5*pi,0.7*pi,0.5,40,3,1)

LVDesignDigFiltViaImpInv(0.5*pi,0.7*pi,0.5,40,4,1)
% Page 182-----------------------------------------------------------------
LVs2zViaBilinearEx(3,1,1,1,1)
LVs2zViaBilinearEx(3,1,1,5,1)
LVs2zViaBilinearEx(3,1,5,5,1)
LVs2zViaBilinearEx(3,1,5,25,1)
% Page 183-----------------------------------------------------------------
LVs2zViaBilinearEx(4,[],4,8,0)
% -------------------------------------------------------------------------
LVs2zViaBilinearEx(3,1,1,5,1)
% Page 184-----------------------------------------------------------------
LVsCheby1LPF2zCheby1LPF(1,40,0.5*pi,0.6*pi)
% Page 186-----------------------------------------------------------------
[b,a] = ellip(7,0.5,50,[0.3,0.5])
Hz = LVzFr(b,a,1000,21);
% Page 187-----------------------------------------------------------------
[b,a] = cheby1(7,0.25,[0.3,0.5],'stop')
Hz = LVzFr(b,a,1500,18);
[rS1,rS2,Rp1,Rp2] = LVxBW4DigitalCheb1Notch(Hz,0.3,0.5,50)
% -------------------------------------------------------------------------
[z,p,k] = ellip(6,0.25,50,[0.33,0.5],'stop')
Hz = LVzFr(k*poly(z),poly(p),1000,17);
% Page 191-----------------------------------------------------------------
Hlp = [2.2,1; 2.05,2.2]; hveclp = [-2.05,-0.883]'
a = pinv(Hlp)*hveclp
a = [1,a']'

Hup = [1,0,0;2.2,1,0];
b = Hup*a
% Page 193-----------------------------------------------------------------
h = [1,2,3,4,5]; b = [1,4,9,16,25]';
H= [h(1)^0, h(1)^1, h(1)^2; h(2)^0, h(2)^1, h(2)^2; ...
h(3)^0, h(3)^1, h(3)^2; h(4)^0, h(4)^1, h(4)^2; ...
h(5)^0,h(5)^1,h(5)^2];
a = b\H
% -------------------------------------------------------------------------
Lim = 50; x = -Lim/2:1:Lim/2; for RowCtr = -Lim/2:1:Lim/2;
XX(RowCtr + Lim/2 + 1,1:1:3) = [1 x(RowCtr + Lim/2 + 1) ...
x(RowCtr+Lim/2+1)^2]; end;
Y = ([-Lim/2:1:Lim/2].^2)'; a = Y\XX
% Page 194-----------------------------------------------------------------
[b, a] = LVxProny([(0.9).^( 0:1:50 )],2,2)





