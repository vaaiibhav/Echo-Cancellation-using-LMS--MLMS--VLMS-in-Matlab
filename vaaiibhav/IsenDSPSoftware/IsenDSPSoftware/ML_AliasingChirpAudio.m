function  ML_AliasingChirpAudio(SR,StartFreq,EndFreq)
% ML_AliasingChirpAudio(3000,0,1500)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

StartFreq = floor(abs(StartFreq));
%if StartFreq<1
%    StartFreq=1;
%end

EndFreq = fix(abs(EndFreq));
if EndFreq<1
    EndFreq=1;
end

ChirpAudioSR = SR;

t = 0:1/ChirpAudioSR:1;
ChirpAudio = chirp(t,StartFreq,1-1/ChirpAudioSR,EndFreq);

figure(87)
subplot(211)
plot(EndFreq*t,ChirpAudio)
xlabel(['A Chirp, ',num2str(StartFreq),' Hz to ',num2str(EndFreq),' Hz, Sampling Rate =  ',...
      num2str(ChirpAudioSR),' Hz; X-Axis = Frequency'])
ylabel(['Amplitude'])
axis([0 EndFreq -inf inf])

ChirpAudio = (0.99/max(abs(ChirpAudio)))*ChirpAudio';

SizeFFT = 64;
Sampsovlap = SizeFFT/2;
inpcntFft = length(ChirpAudio);
%==================================================================
y = zeros(SizeFFT,ceil(inpcntFft/SizeFFT));
y(1:SizeFFT,1) = ChirpAudio(1:SizeFFT,1);
remainder = inpcntFft-(2*SizeFFT-Sampsovlap);
n = 2;
while remainder >= SizeFFT-Sampsovlap;
y(1:SizeFFT,n) = ...
ChirpAudio(((n-1)*SizeFFT+1-Sampsovlap*(n-1)):(n*SizeFFT-Sampsovlap*(n-1)),1);
n=n+1;
remainder = inpcntFft-(n*SizeFFT-(n-1)*Sampsovlap);
end
IntNoWin = n;
y(1:SizeFFT,IntNoWin+1) = zeros(SizeFFT,1);
y(1:remainder,IntNoWin+1) = ...
ChirpAudio(((IntNoWin*SizeFFT)-(n-1)*Sampsovlap+1):inpcntFft,1);

%==================================================================
fry = fft(y);
fry = abs(fry(1:Sampsovlap+1,:));

szfry = size(fry);
maxs = max(fry);
normFreq = zeros(1,szfry(2));
for ctr = 1:1:szfry(2)
    testsig = fry(:,ctr);
    loc = find(testsig==maxs(ctr));
    normFreq(1,ctr) = (loc(1)-1)/Sampsovlap;
end

subplot(212)
% plot([0:1:szfry(2)-1]/szfry(2),normFreq)
stem([0:1:szfry(2)-1]/szfry(2),normFreq)
xlabel('Time, sec')
ylabel('Norm. Freq.')
axis([0, 1, 0, 1])

sound(ChirpAudio, SR)