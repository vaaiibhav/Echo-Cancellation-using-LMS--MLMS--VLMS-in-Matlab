function LVLPFViaSymm2AsymmIDFT(Ak) 
% LVLPFViaSymm2AsymmIDFT([0,0,1,1,1,1,1,0,0]); % odd
% LVLPFViaSymm2AsymmIDFT([0,0,1,1,1,1,1,0,0,0]) % even
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
L = length(Ak); if ~(rem(L,2)==0) % odd length filter
M = (L-1)/2; symmDFT = Ak.*exp(-j*pi*[-M:1:M]*(2*M)/L);
LenNegBins = M; NegBins = symmDFT(1,1:LenNegBins);
ZeroPosBins = symmDFT(1,LenNegBins+1:length(symmDFT));
NetDFT = [ZeroPosBins  NegBins]; Imp = real(ifft(NetDFT));
k = 0:1:(L-1)/2; 
else % even length filter
symmDFT = Ak.*exp(-j*pi*[-L/2+1:1:L/2]*(L-1)/L);
LenNegBins = L/2-1; NegBins = symmDFT(1,1:LenNegBins);   
ZeroPosBins = symmDFT(1,LenNegBins+1:length(symmDFT));
NetDFT = [ZeroPosBins  NegBins]; Imp = real(ifft(NetDFT));
k = 0:1:L/2; 
end
figure(3); clf; subplot(211); stem([0:1:length(Imp)-1],Imp); 
xlabel('Sample'); ylabel('Amplitude')
fr = abs(fft(Imp,1024));subplot(212); 
plot([0:1:512]/512,fr(1,1:513)); 
xlabel('Frequency, Units of \pi')
ylabel('Magnitude'); hold on; for ctr = k;
plot([(2*ctr/L),(2*ctr/L)],[0 1],'b:'); 
plot([2*ctr/L],abs(ZeroPosBins(1,ctr+1)),'ko');
axis([0 1 -inf inf]); end