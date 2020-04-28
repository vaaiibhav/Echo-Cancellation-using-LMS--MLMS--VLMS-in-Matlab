function LVxFRviaCirConDFTs(LL,LenWin,wintype)
% LVxFRviaCirConDFTs(4096,20,1)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
IdealLPF = LVIdealLPFImpResp(0.5*pi,LL);
dft1 = fft(IdealLPF);

figure(66)

if ~(rem(LenWin,2)==0)
LenWin = LenWin + 1;
end
offset = LL/2 - LenWin/2;

thewin = zeros(1,LL);

if wintype==1 % rect
thewin(1,offset+1:LL/2+LenWin/2) = rectwin(LenWin)';
plotLoLim = -50;
elseif wintype==2
thewin(1,offset+1:LL/2+LenWin/2) = hamming(LenWin)';
plotLoLim = -70;
elseif wintype==3
thewin(1,offset+1:LL/2+LenWin/2) = blackman(LenWin)';
plotLoLim = -90;
elseif wintype==4
thewin(1,offset+1:LL/2+LenWin/2) = kaiser(LenWin,10)';
plotLoLim = -110;
end

subplot(311)
xplot = [0:1:LL/2]/(LL/2);

ploty = abs(dft1(1,1:LL/2+1));
ploty = ploty/max(ploty);
plot(xplot,20*log10(ploty+eps))
xlabel('(a) Normalized Frequency')
ylabel('Magnitude')
axis([0 1 0.7*plotLoLim 5])

subplot(312)

dft2 = fft(thewin);
ploty = abs(dft2(1,1:LL/2+1));
ploty = ploty/max(ploty);
plot(xplot,20*log10(ploty+eps))
xlabel('(b) Normalized Frequency')
ylabel('Magnitude')
axis([0 1 plotLoLim 5])

subplot(313)
theCircCon = FastCirCon(dft1,dft2);
finalans = real(ifft(theCircCon));
finalans = finalans/max(abs(finalans));

ploty = abs(theCircCon(1,1:LL/2+1));
ploty = ploty/max(ploty);
plot(xplot,20*log10(ploty+eps))
xlabel('(c) Normalized Frequency')
ylabel('Magnitude')
axis([0 1 plotLoLim 5])

return
figure(9); 

wf = IdealLPF.*thewin;
ansy = find(abs(wf)>0);
wf = wf(ansy);

subplot(211)
stem(wf)
subplot(212)
yyy = abs(fft(wf,2^15));
yyy = yyy(1,1:length(yyy)/2);
yyy = yyy/max(yyy);
plot([0:1:length(yyy)-1]/length(yyy), 20*log10(yyy+eps));
axis([0 1 plotLoLim 5])

figure(7);
if rem(LenWin,2)==0
plotbounds = LL/2-LenWin : LL/2+LenWin;
xxxplot = [0:1:length(plotbounds)-1]-(LenWin/2+1);
hndPlot = stem(xxxplot,finalans(1,plotbounds));
xlabel('Sample')
ylabel('Amplitude')
axis([min(xxxplot)  max(xxxplot) 1.5*min(finalans) 1.1*max(finalans)])
end

function b = LVIdealLPFImpResp(wc,L)
M = (L-1)/2; n = 0:1:L-1;
b = sin(wc*(n - M + eps))./(pi*(n - M + eps));

function cc = FastCirCon(s1,s2)
cc = fft(   real(ifft(s1)).*real(ifft(s2)) );


