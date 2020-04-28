% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
% this script is not meant to be called; it contains clean copies of m-code
% from the text. To run any of the example m-code segments below, copy the
% code, paste it into the Command line, and press Return.
return
%Page 102-----------------------------------------------------------------
x = [(2+j),-1,j,3];N= length(x);
n = 0:1:N-1; W= exp(-j*2*pi/N);
for k = 0:1:N-1, DFS(k +1) = sum(x.*(W.^(n*k))); end
%--------------------------------------------------------------------------
N= length(DFS); k = 0:1:N-1;W= exp(j*2*pi/N);
for n = 0:1:N-1, x(n +1) = sum(DFS.*(W.^(n*k))); end
x = x/N
%--------------------------------------------------------------------------
x = [(2+j),-1,j,3];N= length(x); n = 0:1:N-1;
k = n; DFS = x*(exp(-j*2*pi*(n'*k)/N))
N= length(DFS); n = 0:1:N-1;
k = n; x = DFS*(exp(j*2*pi*(n'*k)/N))/N
%Page 103-------------------------------------------------------------------
x = [1,0,1];N= length(x); n = 0:1:N-1;W= exp(-j*2*pi/N);
k = 0:1:N-1; DFS = x*(W.^(n'*k)),
figure(9); stem(n,abs(DFS))
xlabel('Normalized Frequency'); ylabel('Magnitude')
%--------------------------------------------------------------------------
x = [1 0 1];N= length(x); n = 0:1:N-1; k = n;
W= exp(-j*2*pi/N); DFS = x*(W.^(n'*k));
w = 0:0.01:2*pi; DTFT = 1+exp(-j*2*w);
figure(8); clf; hold on; plot(w/(pi),abs(DTFT));
stem(2*[0:1:2]/3,abs(DFS),'ro');
xlabel('NormFreq (Units of pi)'); ylabel('Magnitude')
%Page 105------------------------------------------------------------------
 LV_TDReconViaSampZXform([ones(1,4)],[1],8,[-15:1:15])
 LV_TDReconViaSampZXform([ones(1,4)],[1],5,[-15:1:15])
 LV_TDReconViaSampZXform([ones(1,4)],[1],3,[-15:1:15])
%--------------------------------------------------------------------------
x = [1,1,1,1]; n = [0,1,2,3];
[y, nOut] = LVAddSeqs(x,n,x,n+3);
[y, nOut] = LVAddSeqs(y,nOut,x,n-3);
[y, nOut] = LVAddSeqs(y,nOut,x,n+6);
[y, nOut] = LVAddSeqs(y,nOut,x,n-6);
figure(7); stem(nOut,y)
%Page 106------------------------------------------------------------------
inc = 1/999; xvec = inc/2:inc:1; zp = 2*pi*xvec; z = exp(j*zp);
ZX1 = 1 + z.^(-1) + z.^(-2) + z.^(-3);
figure(33); subplot(211);
plot(xvec,abs(ZX1)); axis([0,1,0,5])
k = 0:1:3; Fndx = fix((k/4)*1000) + 1;
Xtil = ZX1(Fndx); S = 0;
for k = 0:1:3
S = S + Xtil(k+1)./(1-(exp(j*2*pi*k/4))*(z.^(-1)));
end
ZX2 = S.*((1-z.^(-4))/4);
subplot(212); plot(xvec,abs(ZX2))
axis([0,1,0,5])
%Page 107------------------------------------------------------------------
LVxZxformFromSamps([1,1,1,1],[1],6,5000,30)
%Page 111------------------------------------------------------------------
s = [ones(1,64)]; y = fft(s);
subplot(2,1,1); stem(abs(y)); subplot(2,1,2); stem(unwrap(angle(y)))
%Page 112------------------------------------------------------------------
n = 0:1:7; x = [0:1:7]; xret = x(mod(-n,8)+1);
subplot(3,2,1);stem(n,x); y = fft(x); yret = fft(xret);
subplot(3,2,2);stem(n,xret); subplot(3,2,3); stem(n,real(y));
subplot(3,2,4);stem(n,real(yret));subplot(3,2,5); stem(n,imag(y));
subplot(3,2,6); stem(n,imag(yret))
%--------------------------------------------------------------------------
x=[1 2 3 4]; xsh1 = [4 1 2 3]; xft = fft(x), xshfft = fft(xsh1),
k=0:1:length(x)-1; prxshfft = xft.*(exp(-j*2*pi*1*k/4))
%Page 113------------------------------------------------------------------
N = 9; x = randn(1,N); td = sum(abs(x).^2),
fd = (1/N)*sum(abs(fft(x)).^2)
%--------------------------------------------------------------------------
n=[0:1:31]; x=[16:-1:-15]; y=fft(x); figure
subplot(2,1,1); stem(n,real(y)); subplot(2,1,2); stem(n,imag(y))
%--------------------------------------------------------------------------
n=[0:1:32]; x=[16:-1:-16]; y=fft(x); figure
subplot(2,1,1); stem(n,real(y)); subplot(2,1,2); stem(n,imag(y))
%Page 115-------------------------------------------------------------------
x=[1:1:8]; subx = x(1,2:length(x)); xe = 0.5*(subx + fliplr(subx));
xo = 0.5*(subx - fliplr(subx)); xeven = [x(1), xe], xodd = [0, xo];
DFTevenpt = fft(xeven), ReDFTx = real(fft(x))
N=length(x); var = sum([(DFTevenpt-ReDFTx)/N].^2)
%Page 116------------------------------------------------------------------
cos(2*pi*1.3*(0:1:3)/4) + sin(2*pi*(0.85)*(0:1:3)/4)
%--------------------------------------------------------------------------
fft([cos(2*pi*1.3*(0:1:3)/4) + sin(2*pi*(0.85)*(0:1:3)/4)])
%Page 117-------------------------------------------------------------------
y = ones(1,32); stem(abs(fft(y)))
x = 0:1:length(y) - 1; stem(x, abs(fft(y)))
%Page 118------------------------------------------------------------------
k = 0; y = sum(exp(-j*2*pi*(0:1:2)*k/3 ).*([1,-1,1]))
%Page 121------------------------------------------------------------------
x = [1 2 3 4]; N=length(x); nkvec = 0:1:N-1;
W= exp(nkvec'*nkvec).^(-j*2*pi/N); dft = W*x',
mfft = fft(x)
%--------------------------------------------------------------------------
x = [(1+j) 2 (3+j) 4]; N=length(x); nkvec = 0:1:N-1;
W= exp(nkvec'*nkvec).^(-j*2*pi/N); dft = W*conj(x'),
mfft = fft(x)
%--------------------------------------------------------------------------
LVDFTCompute(0)
%Page 124------------------------------------------------------------------
LVxDFTComputeSawtooth
%Page 125------------------------------------------------------------------
LVxDFTComputeSymmIndex
%Page 126------------------------------------------------------------------
LVxDFTComputeImpUnBal
%Page 127------------------------------------------------------------------
SR = 32; n = 0:1:SR-1; w = 0;
for ctr = 1:1:SR; w = w + (1/SR)*cos(2*pi*n/SR*(ctr-1)); end;
dftans = fft(w); figure(9); subplot(311); stem(w);
subplot(312); stem(real(dftans));
subplot(313); stem(imag(dftans))
%Page 137-------------------------------------------------------------------
y = fft([1, 2, 3, 4])
%--------------------------------------------------------------------------
Signal = [0:1:7];
J = 0;
for I = 1:1:length(Signal)-2
k = length(Signal)/2;
while (k <= J)
J = J - k;
k = k/2;
end
J = J + k;
if I < J
temp = Signal( J+1);
Signal( J+1) = Signal(I+1);
Signal(I+1) = temp;
end
end
BitReversedSignal = Signal
%page 138------------------------------------------------------------------
LV_FFT(8,0)
%--------------------------------------------------------------------------
x = BitReversedSignal;
Rx = real(x); Ix = imag(x);
M = log2(length(x)); LenSig = length(x);
for L = 1:1:M
LE = 2^L; LE2 = LE/2;
uR = 1; uI = 0;
sR = cos(pi/LE2); sI = -sin(pi/LE2);
for J = 1:1:LE2
Jmin1 = J-1;
for I = Jmin1:LE:LenSig - 1
Ip = I + LE2;
tR = Rx(Ip+1)*uR - Ix(Ip+1)*uI;
tI = Rx(Ip+1)*uI + Ix(Ip+1)*uR;
Rx(Ip+1) = Rx(I+1) - tR;
Ix(Ip+1) = Ix(I+1) - tI;
Rx(I+1) = Rx(I+1) + tR;
Ix(I+1) = Ix(I+1) + tI;
end
tR = uR;
uR = tR*sR - uI*sI;
uI = tR*sI + uI*sR;
end
end
x = Rx + j*Ix
%Page 140------------------------------------------------------------------
theBin = 2; x = randn(1,8);N= length(x);
fBin = exp(-j*2*pi*(0:1:N-1)*theBin/N);
cc = fliplr(fBin); y = conv(cc,x);
CZL = y(1,N), mftx = fft(x); fftBin = mftx(theBin+1)
%Page 141------------------------------------------------------------------
p = exp(j*2*pi*1/4); ans1 = p.^(4:-1:1), ans2 = p.^(0:-1:-3)
%--------------------------------------------------------------------------
Bin = 3; x = [-3:1:4]; N=length(x);
GBin = filter(1,[1 -(exp(j*2*pi*Bin/N))],[x 0]);
GoertzelBin = GBin(N+1), ft = fft(x); ftBin = ft(Bin+1)
%Page 142------------------------------------------------------------------
Bin = 3; x = [-3:1:4]; N = length(x);
GB = filter(1,[1 -2*(cos(2*pi*Bin/N)) 1],[x, 0]);
GrtzlBin = GB(N+1) - exp(-j*2*pi*Bin/N)*GB(N),
ft = fft(x); ftBin = ft(Bin+1)
%--------------------------------------------------------------------------
MagSqGBin = GB(N+1)^2 - 2*cos(2*pi*Bin/N)*GB(N+1)*GB(N) + GB(N)^2,
AbsGBin = MagSqGBin^0.5
%Page 145------------------------------------------------------------------
y = mod(-2,8)
%Page 147------------------------------------------------------------------
b = [1,2,2,1]; x = [1,0,-1,1]; N = 4; n = 1:1:N;
for m = 1:1:N; CirCon(m) = sum(b(n).*(x(mod(m-n,N)+1)));
end; subplot(2,1,1), stem(CirCon)
subplot(2,1,2); y = real(ifft(fft(b).*fft(x))); stem(y)
%--------------------------------------------------------------------------
figure(120); subplot(3,1,1); stem(conv( [1,2,2,1], [1,0,-1,1]))
%Page 148------------------------------------------------------------------
b = [1,2,2,1,zeros(1,3)]; x = [1,0,-1,1,zeros(1,3)];
N= 7; n = 1:1:N;subplot(3,1,2)
for m = 1:1:N; CirCon(m) = ...
sum(b(n).*(x(mod(m-n,N)+1))); end;
stem(CirCon); subplot(3,1,3);
y = real(ifft(fft(b,8).*fft(x,8))); stem(y(1:7))
%Page 149------------------------------------------------------------------
SR = 1000; ts = chirp( [0:1/(SR-1):1],0,1,SR/2);
s = [1,1,1,1]; lincon = conv(s, ts);
linconbyfft = real(ifft( fft(ts, 1024).*fft(s,1024)));
figure(14); subplot(2,1,1); plot(lincon);
subplot(2,1,2); plot(linconbyfft(1:1003))
%Page 151------------------------------------------------------------------
FirstConv = real(ifft(fft([1,1,0,0]).*fft([1,2,0,0])))
%Page 152------------------------------------------------------------------
SecConv = real(ifft(fft([1,1,0,0]).*fft([3,4,0,0])))
%--------------------------------------------------------------------------
y = conv([1,1], [1,2,3,4])
%Page 159------------------------------------------------------------------
LVDTFTWindowsOnly(64)
LVDTFTWindowing(64)
%Page 162------------------------------------------------------------------
LVxWindowingDisplay(1,66,68)
%--------------------------------------------------------------------------
LVxWindowingDisplay(1,66.5,68.6)
%Page 163------------------------------------------------------------------
LVxWindowingDisplay(1,66.5,70.7)
%Page 164------------------------------------------------------------------
x = [1,0,0,0,0,0,0,1]; y = fft(x,1024); plot(abs(y))
%Page 165------------------------------------------------------------------
LVDTFTUsingPaddedFFT([1,0,0,0,0,0,0,1],1024)
%Page 166------------------------------------------------------------------
LVPaddedDFTMovie([1,0,0,0,0,0,0,1],128,1)
%Page 169------------------------------------------------------------------
y = abs(fft(conv(conv((0.9*j).^(0:1:99),(-0.9*j).^(0:1:99)),...
[1,0,0,0,0,0,1]),1024)); figure(8); subplot(211);
plot(y(1,1:513)); axis([0,513,0,inf ])
a = conv([1,0.9*j],[1,-0.9*j]); b = [1,0,0,0,0,0,1]
ans = filter(b,a,[1,zeros(1,1000)]);fr=abs(fft(ans,1024));
subplot(212); plot(fr(1,1:513))
axis([0,513,0,inf ])
%Page 171-------------------------------------------------------------------
s = [2,1,-1,1]; n = 0:1:3; F = fft(s); k = [0:1:3]'; Arg = 2*pi*k/4;
for n = 0:1:3; hold on;
stem(n, real(0.25*sum(F(k+1)*exp(j*Arg*n)))); end
%Page 172------------------------------------------------------------------
s = [2,1,-1,1]; n = 0:1:3; F = fft(s); ID = zeros(1,4);
for k=1:1:4; ID = ID+ 0.25*F(k)*exp(j*2*pi*(k-1)*n/4),
pause(1); end
%--------------------------------------------------------------------------
x = [(1+2*j),2,(3-2*j),4];N = length(x); n = 0:1:N-1; k = 0:1:N-1;
MSfft = fft(x);CW = exp(n'*k).^(j*2*pi/N);
idft = (1/N)*CW*conj(MSfft'),MSifft = ifft(MSfft)
%Page 173-------------------------------------------------------------------
idft = (1/4)*conj(fft(conj(fft([1+4*j,2+3*j,3+2*j,4+j]))))
%Page 174------------------------------------------------------------------
LVxInvDFTComputeRect25
%Page 178------------------------------------------------------------------
LVxIDFTviaPosK([randn(1,19)])
%Page 180------------------------------------------------------------------
LVxPhaseShiftViaDFT
%Page 182------------------------------------------------------------------
LVxDFTEqualization(3,0.9,0.2,2048,450)
