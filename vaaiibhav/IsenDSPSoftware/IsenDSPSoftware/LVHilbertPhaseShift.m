function LVHilbertPhaseShift(DestWF,UserWF,SR)
% function LVHilbertPhaseShift(DestWF,UserWF,SR)
%
% DestWF is passed as 1 for a final waveform of sawtooth, 2 for square,
% or 3 to enter a user-supplied waveform as UserWF
% If DestWF is passed as 3, SR may be passed as [], if DestWF is passed
% as 1 or 2, UserWF may be passed as [].
% SR is the desired length of the impulse response to be computed when DestWF
% is passed as 1 or 2.
% Computes Hilbert transform using IDFT with two methods, real FD Hilbert
% mask and imaginary FD Hilber mask and displays the real and imaginary
% parts of the inverse DFT of the product of the DFT of the signal and the
% respective FD Hilbert mask.
%
% Test Calls: 
%
% LVHilbertPhaseShift(1,[],64) % Sawtooth 
% LVHilbertPhaseShift(2,[],64) % Square
% LVHilbertPhaseShift(3,[cos(2*pi*16*[0:1:63]/64)],[])
%  
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool


if DestWF==1|DestWF==2
    WF = zeros(1,SR);
t = (0:1/SR:1-1/SR);
freq = 2;
ctr = 1:DestWF:SR/(2*freq);
for theCtr = ctr
   WF = WF + (1/theCtr)*cos(2*pi*theCtr*freq*t);
end
else % user-supplied test signal
    WF = UserWF;
end
[PhaseShiftedOutputRM] = HilbertPhaseShift(WF,1);
[PhaseShiftedOutputIM] = HilbertPhaseShift(WF,2);

figure(467)
clf
plotlimHi = max([1.2,1.2*max(abs(WF))  ]   );
plotlimLo = min([-1.2,-1.2*max(abs(WF))   ]  );

subplot(321)
xvec = 0:1:length(WF)-1;
if SR>64
plot(xvec,real(WF))
else
stem(xvec,real(WF),'bo');
end
ylabel(['Real'])
xlabel(['(a) Sample, Test Signal'])
axis([0, length(WF), plotlimLo, plotlimHi])

subplot(322)
if SR>64
plot(xvec,imag(WF))
else
stem(xvec,imag(WF),'bo');
end
ylabel(['Imag'])
xlabel(['(b) Sample, Test Signal'])
axis([0, length(WF), plotlimLo, plotlimHi])

subplot(323)
if SR>64
plot(xvec,real(PhaseShiftedOutputIM))
else
stem(xvec,real(PhaseShiftedOutputIM),'bo');
end
xlabel(['(c) Sample, Hilbert via Imag FD Filt'])
ylabel(['Real'])
axis([0, length(WF), plotlimLo, plotlimHi])

subplot(324)
if SR>64
plot(xvec,imag(PhaseShiftedOutputIM))
else
stem(xvec,imag(PhaseShiftedOutputIM),'bo');
end
xlabel(['(d) Sample, Hilbert via Imag FD Filt'])
ylabel(['Imag'])
axis([0, length(WF), plotlimLo, plotlimHi])

% right real FD mask plots

subplot(325)
if SR>64
plot(xvec,real(PhaseShiftedOutputRM))
else
stem(xvec,real(PhaseShiftedOutputRM));
end
xlabel(['(e) Sample, Hilbert via Real FD Filt'])
ylabel(['Real'])
axis([0, length(WF), plotlimLo, plotlimHi])

subplot(326)
if SR>64
plot(xvec,imag(PhaseShiftedOutputRM))
else
stem(xvec,imag(PhaseShiftedOutputRM),'bo');
end
xlabel(['(f) Sample, Hilbert via Real FD Filt'])
ylabel(['Imag'])
axis([0, length(WF), plotlimLo, plotlimHi])


function [PhaseShiftedOutput] = HilbertPhaseShift(WF,FDMethod)

LL = length(WF);  % make FD Mask and TDHilbert for circ con of same length as test signal
sigFD = fft(WF);

if FDMethod==1  % All-Real FD mask
    if rem(LL,2)==0  % even length
   FreqDomHilbert = [ 1 2*ones(1,(LL/2)-1) 1 zeros(1,(LL/2)-1) ];
    else   % odd length
   FreqDomHilbert = [ 1 2*ones(1,((LL-1)/2)) zeros(1,((LL-1)/2)) ];
    end     
 PhaseShiftedOutput = (ifft(sigFD.*FreqDomHilbert)); 
else  % All-Imag FD Mask
    if rem(LL,2)==0
   FreqDomHilbert = zeros(1,LL) + j*([0 -ones(1,(LL/2)-1) 0 ones(1,(LL/2)-1)]);
    else   
   FreqDomHilbert = zeros(1,LL) + j*([0 -ones(1,((LL-1)/2)) ones(1,((LL-1)/2)) ]);
    end
   PhaseShiftedOutput = (ifft(sigFD.*FreqDomHilbert));
end

return