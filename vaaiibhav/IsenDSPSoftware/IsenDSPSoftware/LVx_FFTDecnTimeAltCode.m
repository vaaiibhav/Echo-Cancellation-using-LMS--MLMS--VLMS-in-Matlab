function LVx_FFTDecnTimeAltCode(x,CmprMode,BitRCd,BtrFCd)
% function LVx_FFTDecnTimeAltCode(x,CmprMode,BitRCd,BtrFCd)
% Pass CmprMode as 0 to compare the time of execution for various user-written m-coded
% FFT's to the built-in fft; Pass CmprMode as 1 to compare the time of execution of a
% direct DFT to the built-in fft.
% When CmprMode is passed as 0, pass BitRCd as 0 to use a standard
% bit-reversal routine (from text)  or as 1 to use direct decimation (user-written). Pass BtrFCd as 0
% to run the butterfly routine from the text or as 1 to run a user-written
% butterfly routine. When CmprMode is passed as 1, you may pass BitRCd and BtrFCd as the empty
% matrix [].
% Example calls:
% LVx_FFTDecnTimeAltCode([0:1:31],1,0,0)
% LVx_FFTDecnTimeAltCode([0:1:31],0,0,0)
% LVx_FFTDecnTimeAltCode([0:1:31],0,1,0)
% LVx_FFTDecnTimeAltCode([0:1:31],0,0,1)
% LVx_FFTDecnTimeAltCode([0:1:31],0,1,1)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

LenSig = length(x);

if LenSig > 2^18
   Comment = 'Signal length greater than 2^18 not allowed due to time required'
   return
elseif (LenSig > 2^12) & (CmprMode==1)
   Comment = 'Signal length greater than 2^12 not allowed when using direct DFT due to time required'
   return
end

N = 1;
while LenSig > 2^N
    N = N + 1;
end
if LenSig == 2^N
else
    NeededLength = 2^N;
    x = [x  zeros(1,NeededLength - LenSig)];  
    LenSig = length(x)
end 

xvec = 0:1:LenSig-1;

if LenSig > 2^12
    shortxvec = 0:1:4095;
else
    shortxvec = xvec;
end

halfxvec = [1:1:LenSig/2];

figure(10)
clf

subplot(221)

if LenSig>64
plot(shortxvec,x(shortxvec+1)) 
else   
stem(xvec,x);
end
ylabel(['Amplitude'])

if length(shortxvec) == length(xvec)
xlabel(['(a) x (Test Signal)'])
else
xlabel(['(a) x(0:1:',num2str(length(shortxvec)-1),')'])
end
axis([0  (length(shortxvec)-1)  min([0   1.1*min(x)])   1.2*max(abs(x))])

tic;
y = abs(fft(x));
ElapsedTimeMathScriptFFT = toc;

strText = ['Elapsed Time for Built-In FFT = ',num2str(ElapsedTimeMathScriptFFT,4),'; '];

subplot(224)
y = abs(fft(x));
if LenSig>64
plot(halfxvec-1,y(halfxvec)); 
else  
stem(halfxvec-1,y(halfxvec));
end

ylabel(['Magnitude'])
xlabel(['(d) Built-In FFT (Bins 0-N/2)'])
axis([0  LenSig/2  0   1.2*max(y) ])

M = log2(length(x));

if CmprMode == 0
% script fft routines==================================================

if BitRCd==0
% standard bit-reversal 

tic

J = 0; 
for I = 1:1:LenSig-2     
    k = LenSig/2;     
   while (k <= J)
        J = J - k;
        k = k/2;
    end
        J = J + k;
    if I < J
        temp = x(J+1);
        x(J+1) = x(I+1);
        x(I+1) = temp;
    end
end 

ElapsedTimeStdBitReversal = toc
%ElapsedTimeStdBitReversal = etime(clock,startTime)
strText = [strText,'Elapsed Time for Std Bit Reversal routine = ',num2str(ElapsedTimeStdBitReversal,4),'; '];

else % recursive (inefficient) decimation-in-time of samples=========================
tic
    
for DecCtr = 0:1:M-2
   NumSegments = 2^(DecCtr);
    for SegCtr = 0:1:NumSegments-1
        LenSeg = length(x)/NumSegments;  
        SegStart = SegCtr*LenSeg;
        SegEnd   = (SegCtr + 1)*LenSeg -1;  
        n = 0:1:LenSeg/2 - 1;  
        segEven = x(1,(SegStart + 2*n)+1);
        segOdd  = x(1,(SegStart + 2*n+1)+1); 
        x(1,SegStart+1:SegEnd+1) = [segEven  segOdd];  
    end
end

ElapsedTimeSampleDecimation = toc
strText = [strText,'Elapsed Time for Decimation-in-time = ',num2str(ElapsedTimeSampleDecimation,4),'; '];
end
%========================================================================
subplot(223)
hold on
if LenSig>64
plot(xvec,x) 
else
stem(xvec,x);
end

ylabel(['Amplitude'])

%axis([0  LenSig  min([0   1.1*min(x)])  1.2*max(abs(x))  ])
if length(shortxvec) == length(xvec)
xlabel(['(c) Time Decimation of x'])
else
xlabel(['(c) Dec. of x (1st ',num2str(length(shortxvec)),' samps)'])
end
axis([0  (length(shortxvec)-1)  min([0   1.1*min(x)])   1.2*max(abs(x))])

% standard textbook butterfly routine
if BtrFCd==0
tic

Rx = real(x);
Ix = imag(x);

for L = 1:1:M
    LE = 2^L;
    LE2 = LE/2;
    uR = 1;
    uI = 0;
    sR = cos(pi/LE2);
    sI = -sin(pi/LE2);
    
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
x = Rx + j*Ix;

ElapsedTimeStdFFTButterfly = toc
strText = [strText,'Elapsed time for Std Butterfly routine = ',num2str(ElapsedTimeStdFFTButterfly,4),'.'];
else
    
tic

% a "direct" and inefficient butterfly routine
for ButterflyCtr = 1:1:M
   NumSegments = 2^(M - ButterflyCtr);
   LenSeg = length(x)/NumSegments;
   
for SegCtr = 0:1:NumSegments-1 
   SegStart = SegCtr*LenSeg;
     k = 0:1:LenSeg/2 - 1;
     TwiddleFac = exp(-j*2*pi*k/LenSeg);      
     Tempsig = x(1,k+SegStart+1);
     s = TwiddleFac.*x(1,k+SegStart + LenSeg/2 + 1);
     x(1,k + SegStart + 1) = Tempsig + s;
     x(1,k + SegStart + LenSeg/2 + 1) = Tempsig - s;       
end
end 

ElapsedTimeAltButterfly = toc
strText = [strText,'Elapsed time for User-Written Butterfly routine = ',num2str(ElapsedTimeAltButterfly,4),'.'];

end

else  % compute direct DFT and plot  
[x,thetime] = DirectFFT(x);
strText = [strText,'Elapsed time for Direct DFT = ',num2str(thetime,4),'.'];

end

strText = strText 

subplot(222)
x = abs(x);
if LenSig>64
plot(halfxvec-1,x(halfxvec));
else
stem(halfxvec-1,x(halfxvec));
end

if CmprMode==0
xlabel(['(b) Script FFT (Bins 0-N/2)'])
else
xlabel(['(b) Script DFT (Bins 0-N/2)'])    
end

ylabel(['Magnitude'])
axis([0,LenSig/2,0,1.2*max(abs(x))])

function [bin,ElapsedtimeDirectDFT] = DirectFFT(x)
tic
N = length(x);
n = 0:1:N-1;
for k = 0:1:N-1
    cmplxcorr = exp(-j*2*pi*n*k/N);
    bin(k+1) = sum(x.*cmplxcorr);
end
ElapsedtimeDirectDFT = toc;