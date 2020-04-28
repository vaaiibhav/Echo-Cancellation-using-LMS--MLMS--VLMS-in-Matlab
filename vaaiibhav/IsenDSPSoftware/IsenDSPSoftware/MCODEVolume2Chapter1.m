% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
% this script is not meant to be called; it contains clean copies of m-code
% from the text. To run any of the example m-code segments below, copy the
% code, paste it into the Command line, and press Return.
return
%Page 7--------------------------------------------------------------------
w = 0:0.01:pi; DTFT = 1./(1-0.9*exp(-j*(w)));
figure(8); plot(w/(pi),abs(DTFT));
xlabel('Normalized Frequency'); ylabel('Magnitude')
%Page 8--------------------------------------------------------------------
w = 0:0.01:pi; DTFT = 1+exp(-j*2*w);
figure(8); plot(w/(pi),abs(DTFT));
xlabel('Normalized Frequency'); ylabel('Magnitude')
%Page 9--------------------------------------------------------------------
LV_DTFT_Basic([1,0,1],300,1)
%Page 10-------------------------------------------------------------------
LVxDTFT([1,0,1],[0:1:2],300,2,1,88)
%Page 12-------------------------------------------------------------------
nN = (0:1:100)/100;
LVxDTFT_MS([cos(2*pi*25*nN)],0,exp(j*2*pi*12.5*nN),200,2,2,1)
%Page 14-------------------------------------------------------------------
N=10^3; dw = 2*pi/N; w = -pi:dw:pi*(1-2/N);
DTFT = 1+exp(-j*2*w);
x0 = (1/(2*pi))*sum(DTFT.*exp(j*w*0)*dw)
x1 = (1/(2*pi))*sum(DTFT.*exp(j*w*1)*dw)
x2 = (1/(2*pi))*sum(DTFT.*exp(j*w*2)*dw)
%Page 15-------------------------------------------------------------------
SR = 100; nN = (0:1:SR)/SR;
LVxDTFT_MS([cos(2*pi*25*nN)],0,exp(j*2*pi*SR*nN),200,2,2,1)
%--------------------------------------------------------------------------
nN = (0:1:100)/100; y = exp(j*2*pi*100*nN)
%Page 16-------------------------------------------------------------------
LVxDTFT_MS([1,0,1],0,exp(j*2*pi*1*(0:1:2)/16),100,2,1,2)
%Page 17-------------------------------------------------------------------
y = conv([1,0,1],[1,0,0,-1])
%Page 18-------------------------------------------------------------------
x = [1:1:9]; [xe,xo,m] = LVEvenOddSymmZero(x,[0:1:8]);
d = LVxDTFT(x,[0:1:8],200,2,2,10);
de = LVxDTFT(xe,m,200,2,2,11);
do = LVxDTFT(xo,m,200,2,2,12);
RMS = sqrt((1/200)*sum((d - (de+do)).^2))
%Page 19-------------------------------------------------------------------
exp(j*2*pi*(0:1:31)*5/32)
%--------------------------------------------------------------------------
w = 5*pi/16; DTFT = 1+exp(-j*2*w);
magH = abs(DTFT)
%--------------------------------------------------------------------------
tdMagResp = max(abs(conv([1,0,1],exp(j*2*pi*(0:1:31)*5/32))))
%Page 20-------------------------------------------------------------------
F = 0.53; N = 128; w = F*pi; dtft = 1+exp(-j*2*w);
M = abs(dtft); theang = angle(dtft);
y1 = M*cos(2*pi*(0:1:N-1)*(N/2*F)/N + theang);
y2 = conv([1 0 1],cos(2*pi*(0:1:N-1)*(N/2*F)/N));
stem(y1,'bo'); hold on; stem(y2,'r*')
%Page 21-------------------------------------------------------------------
N=10^3; dw=2*pi/N; w = 0:dw:2*pi-dw;
H = (1+exp(-j*2*w))./(1-0.9*exp(-j*w));
figure; subplot(2,1,1); plot(abs(H)); subplot(2,1,2); plot(angle(H))

