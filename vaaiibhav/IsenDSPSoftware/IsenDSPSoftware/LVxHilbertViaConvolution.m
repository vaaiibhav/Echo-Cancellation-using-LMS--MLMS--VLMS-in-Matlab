function LVxHilbertViaConvolution(TestSeqLen,TestWaveType,FilterLen,TestSigFreq,FDMethod,UseWin)
%
%TestSeqLength is the desired test sequence length; 
%TestWaveType is passed as '1' for sawtooth, and '2' for squarewave; 
%FilterLen is the desired length of the Hilbert transformer impulse response to be constructed directly using the time domain formula;
%TestSigFreq is the fundamental frequency of the (truncated) sawtooth or square
%wave which is computed and used as the test signal;
%FDMethod is passed as 1 to use an All-Real FD Mask, or 2 for
%All-Imaginary; the mask is the same length as the test signal, and is converted using the ifft 
%to a TD Hilbert transformer which is convolved with the test signal.
%Pass UseWin as 0 to use the raw TD Hilbert impulse response, or 1 to
%window it with a Kaiser window, Beta = 5
%
% Call:
% LVxHilbertViaConvolution(64, 1, 19, 5, 2,0)

% TestWaveType : 1 for sawtooth, 2 for squarewave
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

SR = TestSeqLen;
WF = zeros(1,SR);
t = 0:1/SR:1-1/SR;
freq = TestSigFreq;
L = FilterLen;

n = 0:1:L-1;
M = (L-1)/2;

arg = n-M;
Midpt = find(arg==0);
if isempty(Midpt)
   TDHilbert = 2*(sin(pi*arg./2).^2)./(pi*arg);
else
nonMidpt = find(~(arg==0));
TDHilbert(nonMidpt) = 2*((sin(pi*arg(nonMidpt)./2)).^2)./(pi*arg(nonMidpt));
TDHilbert(Midpt) = 0;
end

if UseWin==1
TDHilbert = kaiser(L,5)'.*TDHilbert;
end

%make up cosine-based WF with harmonics of square wave==================================
for ctr = 1:TestWaveType:fix(SR/(2*freq))
   WF = WF + (1/ctr)*cos(2*pi*freq*ctr*t);
end
%=======================================================================================

LL = length(WF);  % make FD Mask and TDHilbert for circ con of same length as test signal

if FDMethod==1  % All-Real FD mask
    if rem(LL,2)==0  % even length
   FreqDomHilbert = [ 1 2*ones(1,(LL/2)-1) 1 zeros(1,(LL/2)-1) ];
    else   % odd length
   FreqDomHilbert = [ 1 2*ones(1,((LL-1)/2)) zeros(1,((LL-1)/2)) ];
    end    
    
else  % All-Imag FD Mask
    if rem(LL,2)==0
   FreqDomHilbert = zeros(1,LL) + j*([0 -ones(1,(LL/2)-1) 0 ones(1,(LL/2)-1)]);
    else   
   FreqDomHilbert = zeros(1,LL) + j*([0 -ones(1,((LL-1)/2)) ones(1,((LL-1)/2)) ]);
    end
end

TDHilb4CirCon = ifft(FreqDomHilbert);

LVCircularConvolution(WF,TDHilb4CirCon)

PhaseShOutput = conv(TDHilbert,WF);

figure(852)
clf

plotlim1 = 1.1*max(abs(WF));

xvec = 0:1:length(WF)-1;
subplot(311)
stem(xvec,WF,'b');
ylabel(['Amplitude'])
xlabel(['(a) Sample, Test Signal'])
axis([0 SR -plotlim1 plotlim1])

plotlim2 = 1.2*max(abs(TDHilbert));

subplot(312)
stem(0:1:length(TDHilbert)-1,TDHilbert,'bo');
ylabel(['Amplitude'])
xlabel(['(b) Sample, Hilbert Transformer Impulse Response'])
axis([0 L -plotlim2 plotlim2])

plotlim3 = 1.1*max(abs(PhaseShOutput));

subplot(313)
stem(0:1:length(PhaseShOutput)-1,PhaseShOutput);
ylabel(['Amplitude'])
xlabel(['(c) Sample, Test Signal Linearly Convolved with Hilbert Transformer'])
axis([0 SR+L -plotlim3 plotlim3])




