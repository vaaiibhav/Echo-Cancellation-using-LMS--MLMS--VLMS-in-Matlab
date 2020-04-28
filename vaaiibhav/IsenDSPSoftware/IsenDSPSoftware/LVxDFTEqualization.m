function LVxDFTEqualization(tstSig,p,k,SR,xplotlim)
% function LVxDFTEqualization(tstSig,p,k,SR,xplotlim)
% Creates a test signal; tstSig=1 gives 32 periods of a
% length-16 square wave followed by 16 samples equal to zero;
% tstSig=2 yields 64 periods of a length-16 square wave;
% tstSig = 3 yields a length-1024 chirp from 0 to 512 Hz.
% p is a single real pole to be used to generate a decaying 
% amplitude profile which waits SR/4 samples of random noise
% of standard deviation k to generate the test channel 
% impulse response. SR is the length of FFT to be used to 
% model the channel impulse response and create the inverse
% filter. 
% xplotlim is the number of samples of each of the relevant signal
% to plot; the relevant signals are the test signal, the distorted
% test signal, and two equalized versions created respectively by
% frequency domain (FD) deconvolution and time domain(TD)
% deconvolution.
% Since the impulse response is random noise with a decaying magnitude,
% every call to the function generates a completely different impulse
% response for the simulated channel, and the distorted test signal
% will have a different appearance.
% LVxDFTEqualization(1,0.9,0.5,2048,250)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

if tstSig==1
   TestSig = ([ones(1,8),-ones(1,8),zeros(1,16)]')*ones(1,32);
   TestSig = TestSig(:)'; 
elseif tstSig==2
   TestSig = ([ones(1,8),-ones(1,8)]')*ones(1,64);
   TestSig = TestSig(:)'; 
elseif tstSig==3 
   TestSig = chirp([0:1/1023:1],0,1,512);
end

lenTestSig = length(TestSig)

imp = filter(1,[1,-p],[1,zeros(1,SR/4)]);
imp = imp.*(k*randn(1,length(imp)));

disTestSig = filter(imp,1,TestSig);
distsRange = abs(max(disTestSig)-min(disTestSig))
tsRange = abs(max(TestSig)-min(TestSig))
ftsz1 = 13;
ftsz2 = 13;

figure(765)
clf
set(765,'color',[1,1,1])
subplot(221)

if xplotlim > length(TestSig)
    Dxplotlim = length(TestSig);
else
    Dxplotlim = xplotlim;
end
    
plot(TestSig(1,1:Dxplotlim),'k')
xlabel('(a) Test Signal')
ylabel('Amp')
axis([0,xplotlim,[min(TestSig)-0.1*tsRange],[max(TestSig)+0.1*tsRange]])

subplot(222)
plot(disTestSig(1,1:Dxplotlim),'k')
xlabel('(b) Distorted Test Signal')
ylabel('Amp')
axis([0,xplotlim,[min(disTestSig)-0.1*distsRange],[max(disTestSig)+0.1*distsRange]])

frImp = fft(imp,SR);

frEqual = 1./frImp;
TDEqualFilt = real(ifft(frEqual));

figure(765)

lenTDEqFilt = length(TDEqualFilt);

lendisTestSig = length(disTestSig);
lim = max([lenTDEqFilt,lendisTestSig]);
szdisTestSig = size(disTestSig);

netLenZeros = max([(lendisTestSig-lenTDEqFilt),0]);

TDEqualFilt = [TDEqualFilt,zeros(1,netLenZeros)];

disTestSig = [disTestSig, zeros(1,lenTDEqFilt-lendisTestSig)];   
TDViaCircCon = CircConvFCN(TDEqualFilt,disTestSig);

frTestSig = fft(disTestSig,SR);
EqualFR = frTestSig.*frEqual;
EqualTD = real(ifft(EqualFR));

figure(765)

subplot(223)
plot(EqualTD(1,1:Dxplotlim),'k')
xlabel('(c) Equal. Test Sig (Via Inverse FD Filt)')
ylabel('Amp')
axis([0,xplotlim,[min(TestSig)-0.2*tsRange],[max(TestSig)+0.2*tsRange]])

subplot(224)
plot(TDViaCircCon(1,1:Dxplotlim),'k')
xlabel('(d) Equal. Test Sig (Via TD Circ Conv)')
ylabel('Amp')
axis([0,xplotlim,[min(TestSig)-0.2*tsRange],[max(TestSig)+0.2*tsRange]])

function theCircCon = CircConvFCN(Seq1,Seq2)
% CircConvFCN([1 1 1 0 0 0 0 0],[110 121 115 132 115 107 143 117])
% CircConvFCN([1 1 1 1 -1 -1 -1 -1 zeros(1,8)],[0.25 0.5 0.75 1 1 0.75 0.5 0.25 zeros(1,8)])
% CircConvFCN([1 1 1 1 -1 -1 -1 -1],[0.25 0.5 0.75 1 1 0.75 0.5 0.25])
% CircConvFCN([1 1 1 0 0 0 0 0],[110+i*105 121 115 132 115 107 i*143 117])
% CircConvFCN([1 1*i 1 0 0 0 0 0],[110+i*105 121 115 132 115 107 i*143 117])
if ~length(Seq1)==length(Seq2)
   a = 'Sequences Not Same Length!'
   return
end

SR = length(Seq1);
Seq2 = [Seq2(1,1:1),fliplr(Seq2(1,2:length(Seq2)))];

MCircConv = zeros(1,length(Seq1));
maxplotSeq1 = 1.2*max(abs(Seq1));
maxplotSeq2 = 1.2*max(abs(Seq2));
maxplotCircConv = 1.2*max(abs(conv(Seq1,Seq2)));
lagctr = (-length(Seq1)/2) + 1;

for ctr = 1:1:length(Seq2)  
MCircConv(1,ctr:ctr) = sum(Seq1.*(Seq2));

xxx = Seq2(1,length(Seq2):length(Seq2)); % rotates to the right, or counterclockwise
Seq2(1,2:length(Seq2)) = Seq2(1,1:length(Seq2)-1);
Seq2(1,1:1) = xxx;

lagctr =lagctr + 1;

end
theCircCon =  MCircConv;



