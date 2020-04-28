function  LVxPhaseShiftViaDFT
%function  LVxPhaseShiftViaDFT
%
% Creates a signal having the harmonic spectrum of a square wave, but with random
% phases, then resets all phases to 0 degrees initially using the DFT/IDFT, resulting in a
% cusped waveform, which is shown in the second subplot; thereafter, the phase of each 
% bin is shifted a small amount via DFT every time any key is pressed. Eventually, all 
% frequencies in the cusped waveform are shited 90 degrees, and the cusped waveform 
% is converted in a square wave. The third subplot is always the Hilbert
% transform of the waveform in the second plot, and thus starts out as a
% square wave and is gradually converted into a cusped waveform as the
% cusped waveform in the second subplot is converted into a square wave.
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

SR = 1024;
WF = zeros(1,SR);
t = 0:1/SR:1-1/SR;
freq = 3;
for ctr = 1:2:SR/(2*freq)
   WF = WF + (1/ctr)*cos(2*pi*ctr*freq*t + randn(1,1));
end
L = length(WF);
y = fft(WF);

NeutY = (y.*conj(y)).^0.5;

figure(13357)
clf

Comment = 'Press any key for next computation'

for degrees = 0:5:90
PhaseShiftAmt = exp(j*2*pi*(degrees/360));
if rem(L,2)==0
 FreqDomPhiShift = [NeutY(1,1:1),NeutY(1,2:L/2+1)*(PhaseShiftAmt),NeutY(1,L/2+2:L)*conj(PhaseShiftAmt)];
else   
 FreqDomPhiShift = ([NeutY(1,1:1),NeutY(1,2:(L-1)/2+1)*(PhaseShiftAmt),NeutY(1,(L-1)/2+2:L)*conj(PhaseShiftAmt)]);
end
TDResult = real(ifft(FreqDomPhiShift));

plotlim1 = 1.2*max(abs(WF));

subplot(311)
plot(WF)
ylabel(['Amplitude'])
xlabel(['(a) Test Waveform: Cosines w/Random Phase, Amplitudes = 1/(Harmonic No.)'])
axis([0 length(WF) -plotlim1 plotlim1])

plotlim2 = 1.2*max(abs(TDResult));

subplot(312)
plot(TDResult)
ylabel(['Amplitude'])
xlabel(['(b) Signal at (a), all Phases Set to Zero Via DFT/IDFT; then Offset Another ',num2str(degrees),' Degrees'])
axis([0 length(WF) -plotlim2 plotlim2])

%[PhaseShiftedOutput] = imag(hilbert(TDResult));

PhaseShiftAmt = exp(j*2*pi*(90/360));
if rem(L,2)==0
 FreqDomPhiShift = [FreqDomPhiShift(1,1:1),FreqDomPhiShift(1,2:L/2+1)*(PhaseShiftAmt),FreqDomPhiShift(1,L/2+2:L)*conj(PhaseShiftAmt)];
else   
 FreqDomPhiShift = ([FreqDomPhiShift(1,1:1),FreqDomPhiShift(1,2:(L-1)/2+1)*(PhaseShiftAmt),FreqDomPhiShift(1,(L-1)/2+2:L)*conj(PhaseShiftAmt)]);
end
[PhaseShiftedOutput] = real(ifft(FreqDomPhiShift));

plotlim3 = 1.2*max(abs(PhaseShiftedOutput));

subplot(313)
plot(PhaseShiftedOutput)
ylabel(['Amplitude'])
xlabel(['(c) Signal at (b) Shifted 90 Deg Via DFT/IDFT; PRESS ANY KEY TO CONTINUE'])
axis([0 length(PhaseShiftedOutput) -plotlim3 plotlim3])

if degrees==90
    break
end

pause

end
