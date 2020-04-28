function LV_FFT(L, DecOrBitReversal)
% function LV_FFT(L, DecOrBitReversal)
%
% Displays the signal 0:1:L-1 and its progressive decimation using either
% direct decimation or bit reversal. After decimation is complete,
% the computation of the dft via butterflies is progressively displayed.
% Pass DecOrBitReversal = 0 for direct decimation of the test signal, or
% pass DecOrBitReversal = 1 to use bit reversal
%
% Typical calls:
%
% LV_FFT(8, 0)
% LV_FFT(8, 1)
% LV_FFT(256, 0)
% LV_FFT(256, 1)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
Signal = [0:1:L-1];
LenSig = length(Signal);

timeint = 2/L;
N = 1;
while LenSig > 2^N
    N = N + 1;
end

if LenSig == 2^N
else
    NeededLength = 2^N;
    Signal = [Signal  zeros(1,NeededLength - LenSig)];  
    LenSig = length(Signal);
end    
xvec = 0:1:LenSig-1;

figure(10)
clf
ftsz1 = 16;

subplot(2,2,1)
stem(xvec,Signal)
ylabel(['Amplitude'])
xlabel(['(a) Signal'])
axis([0  LenSig  min([0   1.1*min(Signal)])   1.2*max(abs(Signal))])

subplot(2,2,2)
y = abs(fft(Signal));
stem(xvec,abs(fft(Signal)))
ylabel(['Magnitude'])
xlabel(['(b) Built-in FFT (Bin)'])
axis([0  LenSig  0   1.2*max(y) ])

subplot(2,2,4)
xlabel(['(d) FFT (Bin)'])
axis([0  LenSig  0   1.2*max(abs(Signal))])

M = log2(length(Signal));

if DecOrBitReversal==0
for DecCtr = 0:1:M-2
   NumSegments = 2^(DecCtr);
    for SegCtr = 0:1:NumSegments-1
        LenSeg = length(Signal)/NumSegments;  
        SegStart = SegCtr*LenSeg;
        SegEnd   = (SegCtr + 1)*LenSeg -1;  
        n = 0:1:LenSeg/2 - 1;  
        segEven = Signal(1,(SegStart + 2*n)+1);
        segOdd  = Signal(1,(SegStart + 2*n+1)+1); 
        Signal(1,SegStart+1:SegEnd+1) = [segEven  segOdd];
        
        pause(timeint)
            subplot(2,2,3)
            stem(xvec,Signal) 
            ylabel(['Amplitude'])
            xlabel(['(c) Time Decimation of Signal'])
            axis([0  LenSig  min([0   1.1*min(Signal)]) 1.2*max(abs(Signal))  ])
            
    end
end

else
J = 0; 
for I = 1:1:LenSig-2     
    k = LenSig/2;     
   while (k <= J)
        J = J - k;
        k = k/2;
    end
        J = J + k;
    if I < J
        temp = Signal(J+1);
        Signal(J+1) = Signal(I+1);
        subplot(2,2,3)
            stem(xvec,Signal)  % decimated in time
            ylabel(['Amplitude'])
            xlabel(['(c) Time Decimation of Signal'])
            axis([0  LenSig  min([0   1.1*min(Signal)])  1.2*max(abs(Signal))  ])
        
            pause(timeint)
            
            Signal(I+1) = temp;
            subplot(2,2,3)
            stem(xvec,Signal)  % decimated in time
            ylabel(['Amplitude'])
            xlabel(['(c) Time Decimation of Signal'])
            axis([0  LenSig  min([0   1.1*min(Signal)])  1.2*max(abs(Signal))  ])
        
            pause(timeint)
    end
end     
end

subplot(2,2,3)
stem(xvec,Signal)  % decimated in time
ylabel(['Amplitude'])
xlabel(['(c) Time Decimation of Signal'])
axis([0  LenSig  min([0   1.1*min(Signal)])  1.2*max(abs(Signal))  ])

for ButterflyCtr = 1:1:M
   NumSegments = 2^(M - ButterflyCtr);
   LenSeg = length(Signal)/NumSegments;
for SegCtr = 0:1:NumSegments-1 
   SegStart = SegCtr*LenSeg;
      k=0:1:LenSeg/2 - 1;
      TwiddleFac = exp(-j*2*pi*k/LenSeg);
 TempSig1 = Signal(1,k+SegStart+1) + TwiddleFac.*Signal(1,k+SegStart + LenSeg/2 + 1);
 TempSig2 = Signal(1,k+SegStart+1) - TwiddleFac.*Signal(1,k+SegStart+LenSeg/2 + 1);
 Signal(1,k + SegStart + 1) = TempSig1;
 Signal(1,k + SegStart + LenSeg/2 + 1) = TempSig2;       
 
  subplot(2,2,4)
    stem(xvec,abs(Signal))
    ylabel(['Magnitude'])
    xlabel(['(d) FFT (Bin)'])
    axis([0  LenSig  0   1.2*max(abs(Signal))])
    pause(timeint)
 end
end

subplot(2,2,4)
stem(xvec,abs(Signal))
ylabel(['Magnitude'])
xlabel(['(d) FFT (Bin)'])
axis([0  LenSig  0   1.2*max(abs(Signal))])