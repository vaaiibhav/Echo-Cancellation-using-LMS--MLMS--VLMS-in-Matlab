function LVxCyclicVLinearConv(ImpRespSR,SR)
% function LVxCyclicVLinearConv(ImpRespSR,SR)
% ImpRespSR sets the length of an impulse response which
% is made up exclusively of ones, which is then windowed with
% a Hamming window to form a smooth lowpass impulse response.
% SR is the length a of chirp which is convolved with the lowpass
% impulse response using the DFT with padded sequences to effect linear 
% convolution. 
% Test calls:
% LVxCyclicVLinearConv(8,80)
% LVxCyclicVLinearConv(32, 1024)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
if ImpRespSR>SR
   SR = ImpRespSR;
elseif ~(SR/ImpRespSR==floor(SR/ImpRespSR))
   SR = ceil(SR/ImpRespSR)*ImpRespSR;
end

% SR = 2^8;
dispszXL = 12;
t = 0:1/SR:1-1/SR;
tImp = 0:1/ImpRespSR:1 - 1/ImpRespSR;

% Seq1 = sin(2*pi*t*Freq);
% Seq1 = [ones(1,12)  zeros(1,11)  -ones(1,7) 0  ones(1,5)];
% Seq2 = cos(2*pi*t*2*Freq);

% Seq1 = [1 0 1 1 1 -1 -1 1];
% Seq2 = [-1 -1 0 0 1 1 0 -1];
TestSeq = chirp(t,0,1-1/SR,SR/2);
% TestSeq = cos(2*pi*t*20);

% ImpResp = ([ones(1,ImpRespSR)]);
% TestSeq = chirp(t,0,1-1/SR,SR/2);
 %ImpResp = cos(2*pi*tImp*ImpRespSR/4);
N = ImpRespSR;
% ImpResp = (hamming(ImpRespSR)').*([ones(1,ImpRespSR)]);
% ImpResp = (hamming(ImpRespSR)').*(-1).^(0:1:ImpRespSR-1);
 ImpResp = (hamming(ImpRespSR)').*cos(2*pi*(N/4)*[0:1:N-1]/N);

CompositeOutput = zeros(1,length(TestSeq) + length(ImpResp) -1);

LinConv = conv(TestSeq,ImpResp);

LenTestSeq = SR;
LenImpResp = length(ImpResp);
LenTestSubSeq = LenImpResp;

RawSumOfLens = (LenTestSeq+LenImpResp);

pltcrit = length(CompositeOutput);

figure(1200)

maxplotLinConv = max(abs(LinConv));
maxplotTestSeq = max(abs(TestSeq));

subplot(321)
if pltcrit<250
   stem(TestSeq,'bo');
else
   plot(TestSeq)
end

xlabel(['(b) Test Sequence']);
ylabel(['Amplitude']);
axis([0 length(TestSeq)+1 -1.2*maxplotTestSeq 1.2*maxplotTestSeq])

lenImpResp = length(ImpResp);

subplot(321)
if lenImpResp<64
  stem(ImpResp,'bo');
else
   plot(ImpResp)
end

maxplotImpResp = max(abs(ImpResp));

xlabel(['(a) Impulse Response']);
ylabel(['Amplitude']);
axis([0 lenImpResp+1 -1.2*maxplotImpResp 1.2*maxplotImpResp])

LenFFTToUse = 2^(ceil(log2(RawSumOfLens)));
%CompositeOutput(1, (length(TestSeq) + LenImpResp -1 )) = zeros;

HowManySubSeq =  length(TestSeq)/LenTestSubSeq;
IntegerHowMany = floor(length(TestSeq)/LenTestSubSeq);

if IntegerHowMany == HowManySubSeq;    
TheEnd = IntegerHowMany;
else
PadNo = (IntegerHowMany+1)*LenTestSubSeq - length(TestSeq);
TestSeq = [TestSeq zeros(1,PadNo)];
TheEnd = floor(length(TestSeq)/LenTestSubSeq);
end

for ctr = 1:1:TheEnd

ThisTestSubSeq = TestSeq(1,(ctr-1)*LenTestSubSeq+1:ctr*LenTestSubSeq);

%-----replot test signal with ThisTestSubSeq plotted with large circles--------------

subplot(322)

Firstxvec = 1:1:(ctr-1)*LenTestSubSeq;
Midxvec = (ctr-1)*LenTestSubSeq+1:1:ctr*LenTestSubSeq;
Endxvec = ctr*LenTestSubSeq+1:1:length(TestSeq);

if pltcrit<250
   stem(Firstxvec,TestSeq(1,1:(ctr-1)*LenTestSubSeq),'bo');
else
   plot(Firstxvec,TestSeq(1,1:(ctr-1)*LenTestSubSeq),'b')
end

hold on
stem(Midxvec,ThisTestSubSeq,'rd');

if pltcrit<250
   stem(Endxvec,TestSeq(1,ctr*LenTestSubSeq+1:length(TestSeq)),'bo');
else
	plot(Endxvec,TestSeq(1,ctr*LenTestSubSeq+1:length(TestSeq)),'b')
end
hold off

xlabel(['(b) Test Seq, Subseq ',num2str(ctr),' (Diamonds)']);
ylabel(['Amplitude']);
axis([0, length(TestSeq)+1, -1.2*maxplotTestSeq, 1.2*maxplotTestSeq])
%------------------------------------------------------------------------------------

FTPadSeq1 = fft(ThisTestSubSeq,LenFFTToUse);
FTPadSeq2 = fft(ImpResp,LenFFTToUse);
Prod = (FTPadSeq1).*(FTPadSeq2);

padThisSubSeq = [ThisTestSubSeq, zeros(1,length(ThisTestSubSeq))];
padImpResp = [ImpResp, zeros(1,length(ImpResp))];
%-----------------------------------------------------------------------
% middle left plot, padded impulse response

lenPadImpResp = length(padImpResp);

subplot(323)
if lenImpResp<64
   stem(padImpResp,'bo');
else
   plot(padImpResp,'b')
end

xlabel(['(c) Padded Impulse Resp']);
ylabel(['Amplitude']);
axis([0, lenPadImpResp+1, -1.2*maxplotImpResp, 1.2*maxplotImpResp])

%------------------------------------------------------------------------
lenPadSubSeq = length(padThisSubSeq);

subplot(324)
if lenImpResp<64
   stem(padThisSubSeq,'bo');
else
   plot(padThisSubSeq,'b')
end

xlabel(['(d) Padded SubSeq ',num2str(ctr)]);
ylabel(['Amplitude']);
axis([0, lenPadSubSeq+1, -1.2*maxplotTestSeq, 1.2*maxplotTestSeq])

%------------------------------------------------------------------------
EquivLinConv = real(ifft(Prod));

ThingToConcatenateToRight = EquivLinConv(1,1:LenTestSubSeq+LenImpResp-1);

subplot(325)
if pltcrit<64
   stem(LinConv,'bo');
else
   plot(LinConv,'b')
end

xlabel(['(e) Equiv Linear Conv']);
ylabel(['Amplitude']);
% axis([0, inf,     -1.2*maxplotLinConv, 1.2*maxplotLinConv])
axis([0, pltcrit, -1.2*maxplotLinConv, 1.2*maxplotLinConv])

subplot(326)

% this plots composite output before being updated for this iteration
if ctr>1
if pltcrit<64
   stem(CompositeOutput(1,1:(ctr-1)*LenTestSubSeq),'bo');
else
   plot(CompositeOutput(1,1:(ctr-1)*LenTestSubSeq),'b')
end
end

xvec = [];
xvec = (ctr-1)*LenTestSubSeq + 1:1:(ctr-1)*LenTestSubSeq + (LenTestSubSeq);
%---new code---------------------------------------------------------------------------
hold on
if pltcrit<64
   stem(xvec,CompositeOutput(1,(ctr-1)*LenTestSubSeq +1:ctr*LenTestSubSeq),'ro');
else
plot(xvec,CompositeOutput(1,(ctr-1)*LenTestSubSeq+1:ctr*LenTestSubSeq),'r')
end
% superimpose the new segment without adding it numerically----------------------------
if pltcrit<64
   stem((ctr-1)*LenTestSubSeq+1:1:(ctr-1)*LenTestSubSeq +(LenTestSubSeq+LenImpResp-1), ThingToConcatenateToRight,'ko');
else
plot((ctr-1)*LenTestSubSeq+1:1:(ctr-1)*LenTestSubSeq +(LenTestSubSeq+LenImpResp-1), ThingToConcatenateToRight,'k:');
end

xlabel(['(f) Lin Conv Via FFT to Subseq ',num2str(ctr)]);
ylabel(['Amplitude']);
axis([0, pltcrit, -1.2*maxplotLinConv, 1.2*maxplotLinConv])
hold off
pause(1)

% now recompute composite output as final value for this iteration
CompositeOutput(1,(((ctr-1)*(LenTestSubSeq))+1):ctr*(LenTestSubSeq)+LenImpResp-1) = CompositeOutput(1,(((ctr-1)*(LenTestSubSeq))+1):ctr*(LenTestSubSeq)+LenImpResp-1) + ThingToConcatenateToRight;

subplot(326)

line([0 pltcrit],[0 0],'color',[0 0 1]);


if pltcrit<64
   stem(CompositeOutput(1,1:(ctr)*LenTestSubSeq),'bo');
else
plot(CompositeOutput(1,1:(ctr)*LenTestSubSeq),'b')
end

xvec = [];
%----------------------------------
hold on
xvec = (ctr)*LenTestSubSeq + 1:1:(ctr)*LenTestSubSeq + (LenImpResp-1);
if pltcrit<64
   stem(xvec,CompositeOutput(1,(ctr)*LenTestSubSeq +1:ctr*LenTestSubSeq + LenImpResp -1),'ro');
else
plot(xvec,CompositeOutput(1,(ctr)*LenTestSubSeq+1:ctr*LenTestSubSeq + LenImpResp -1),'r')
end
%----------------------------------
xlabel(['(f) Lin Conv Via FFT to Subseq ',num2str(ctr)]);
ylabel(['Amplitude']);
axis([0, pltcrit, -1.2*maxplotLinConv, 1.2*maxplotLinConv])
hold off
pause(1)

end






