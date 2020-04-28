function WF = LVxFIRViaWholeSines(AkPos,AkLOver2,L)
% function WF = LVxFIRViaWholeSines(AkPos,AkLOver2,L) 
%
% Returns an impulse response WF for a Type III or Type IV linear phase FIR.
% AkPos are sample amplitudes for frequencies 1 to (L-1)/2 for odd
% length, or 1 to L/2-1 for even length. AkLOver2 is passed as 0 or the empty matrix []
% for odd values of L, and as a desired amplitude for even values of
% L, which is is the desired filter length
%
% Test calls:
% WF = LVxFIRViaWholeSines([ones(1,38)],[],77); % Hilbert, Type III
% WF = LVxFIRViaWholeSines([0.8,ones(1,36),0.55],[],77); transition values
% aproximately optimized to reduce ripple.
%
% WF = LVxFIRViaWholeSines([1:1:11]*pi/12,[12]*pi/12,24); % differentiator, Type IV
% WF = LVxFIRViaWholeSines([1:1:38]*pi/39,[39]*pi/39,78); % differentiator, Type IV
% WF = LVxFIRViaWholeSines([1:1:38]*pi/39,[],77); % differentiator, Type III
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

%-------------------------------------------------------------------------------------
LenFFT = 1024; M = (L-1)/2;
if rem(L,2)==0
Type = 4; Ak = [AkPos,AkLOver2]; pLenAk = L/2; limK = L/2;
else
Type=3; Ak = [AkPos]; pLenAk = M; limK = M;
end
if length(Ak) < pLenAk
   error('Length of Ak too short for stated value of L')
elseif length(Ak) > pLenAk
   error('Length of Ak too long for stated value of L')
end 
WF = zeros(1,L); n = 0:1:L-1;
for k = 1:1:limK
  if Type==4 & k==limK
       C = 1;
  else
       C = 2;
  end
    newComp = C*Ak(k)*sin(2*pi*(n-M)*k/L);
    WF = WF + newComp;
end   
WF = WF/L; fr = fft(WF,LenFFT);
fr = fr(1,1:LenFFT/2+1)/max(abs(fr));

figure(1997)
clf
subplot(221)
stem(0:1:length(WF)-1,WF,'bo');
plotlim4 = max(abs(WF));
grid on; xlabel(['(a) Sample'])
ylabel(['Amp'])
axis([-1 length(WF) -1.5*abs(min(WF)) 1.2*plotlim4])
bottom = 10^(-5);

subplot(222)
xvec = (0:1:LenFFT/2)/(LenFFT/2);
plot( xvec,abs(fr))
hold on; k = 0:1:limK;
EquivDTFTFreqs = (2*k/L);
theDFT = abs(fft(WF));
theDFT = theDFT/(max(theDFT));
plot(EquivDTFTFreqs, theDFT(k+1),'bo')
hold off; grid on
ylabel(['Mag']); xlabel(['(b)  Frequency'])
axis([0 1 0 1.2])

subplot(223)
plot(xvec,unwrap(angle(fr)))
grid on; xlabel(['(c) Frequency'])
ylabel(['Radians']); axis([0 1  -inf inf])

subplot(224)
xx = 20*log10(abs(fr)+bottom);
plot(xvec,xx)
grid on; xlabel(['(d)  Frequency'])
ylabel(['Mag, dB']); axis([0 1  -100 10])
%==========================================================================


      