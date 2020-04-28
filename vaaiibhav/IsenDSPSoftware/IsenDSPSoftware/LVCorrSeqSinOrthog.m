function LVCorrSeqSinOrthog(LenSeq1,LenSeq2,Freq1,Freq2,phi)
% function DemoCorrSeqSinOrthog(LenSeq1,LenSeq2,Freq1,Freq2,phi)
% LenSeq1 is the number of samples in the first sequence.
% LenSeq2 is the number of samples in the second seq.
% Freq1 is the frequency of the sinusoid in the first seq.
% Freq2 is the frequency of the sinusoid in the second seq.
% phi is a phase offset in radians for the argument of the sinusoid of seq2.
% Test calls:
% LVCorrSeqSinOrthog(128,128,8,8,0)
% LVCorrSeqSinOrthog(32,128,1,8,0) 
% LVCorrSeqSinOrthog(11,55,1,10,0) 
% 
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
LenSeq1 = 11;
LenSeq2 = 55;
Freq1 = 1;
Freq2 = 10;
phi = 0;
if LenSeq1 < 2
    LenSeq1 = 2;
end

if LenSeq2 < 2
    LenSeq2 = 2;
end

if LenSeq2<LenSeq1
    tmpSeq = LenSeq2;
    LenSeq2 = LenSeq1;
    LenSeq1 = tmpSeq;
end

t1 = 0:1/LenSeq1:1-1/LenSeq1;
t2 = 0:1/LenSeq2:1-1/LenSeq2;

Seq1 = sin(2*pi*t1*Freq1);
Seq2 = sin(2*pi*t2*Freq2 + phi);

figure(6018)

subplot(311)
stem(Seq1,'bo');
xlabel(['(a)  First Sequence (Impulse Response); X-Axis = Sample'])
ylabel(['Amplitude'])
axis([0,(LenSeq1+LenSeq2-1),-1.5,1.5])

subplot(312)
stem(Seq2,'bo');
xlabel(['(b)  Second Sequence (Signal); X-Axis = Sample'])
ylabel(['Amplitude'])
axis([0,(LenSeq1+LenSeq2-1),-1.5,1.5])

corrOut = xcorr(Seq2,Seq1);
limmag = max([max(abs(corrOut)),1]);

subplot(313)

startX = length(Seq2)-length(Seq1)+1;
endX = length(corrOut);
lenNetSigOut = endX-startX;
plotlimX = max([(lenNetSigOut+1),10]);

stem(corrOut(1,startX:endX),'bo');
xlabel(['(c)  Correlation Sequence; X-Axis = Sample'])
ylabel(['Amplitude'])
axis([0,plotlimX,-1.2*limmag,1.2*limmag])

%======================================================================
