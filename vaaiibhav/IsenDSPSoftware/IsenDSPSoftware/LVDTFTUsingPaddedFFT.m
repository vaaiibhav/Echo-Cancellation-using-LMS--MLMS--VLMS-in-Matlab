function LVDTFTUsingPaddedFFT(FIRXferFcn,FFTLength)
% LVDTFTUsingPaddedFFT(FIRXferFcn,FFTLength)
%
% FIRXferFcn is a desired impulse response
% which is padded with zeros by the script to the length of the FFT
% which is input as FFTLength.  FFT length is generally specified as many
% times the length of the impulse response.
%
% Typical calls:
%
% LVDTFTUsingPaddedFFT(ones(1,32),1024)
% LVDTFTUsingPaddedFFT(hamming(32)'.*ones(1,32),1024)
% LVDTFTUsingPaddedFFT([1 1 1 1],1024)
% LVDTFTUsingPaddedFFT([1 0 0 1],1024)
% LVDTFTUsingPaddedFFT([1 0 0 0 0 0 0 1],1024)
% LVDTFTUsingPaddedFFT(hamming(33)'.*cos(2*pi*2*(0:1/32:1))+hamming(33)'.*cos(2*pi*8*(0:1/32:1)),1024)
% LVDTFTUsingPaddedFFT(hamming(33)'.*cos(2*pi*7*(0:1/32:1))+hamming(33)'.*cos(2*pi*8*(0:1/32:1)),1024)
% LVDTFTUsingPaddedFFT(exp(-(0:0.2:5.8)),1024)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool


TDSig = zeros(1,60);
TDSig(1,20:29) = ones(1,10);
TDSig(1,30:39) = -ones(1,10);

x = FIRXferFcn;
y = abs(fft(x));
Ny = length(x);
yabsc = 0:1:Ny/2;
yGood = abs(fft(x,FFTLength));
N = FFTLength;
yGoodabsc = 0:1:N/2;

figure(2777)
clf

subplot(311)
stem(0:1:length(x)-1,x,'bo');
xlabel(['(a)  FIR Impulse Response; X-Axis = Sample'])
ylabel(['Amplitude'])
axis([0 length(x) -1.2*abs(max(x)) 1.2*abs(max(x))])

subplot(312)
stem(0:1:length(y)-1,y,'bo');
xlabel(['(b)  ',num2str(Ny),'-pt DFT; X-Axis = DFT Bin'])
ylabel(['Magnitude'])
axis([0 length(x) -inf 1.2*abs(max(y))])

subplot(313)

plot(yGood)
hold on
plot((0:1:length(x)-1)*FFTLength/length(x),y,'bo');
plot((0:1:length(x)-1)*FFTLength/length(x),y,'b:');
xlabel(['(c)  ',num2str(N),'-pt DTFT, Solid; ',num2str(Ny),'-pt DFT, Stars; X-Axis = DTFT Bin'])
ylabel(['Magnitude'])
axis([0 inf -inf 1.2*abs(max(yGood))])



