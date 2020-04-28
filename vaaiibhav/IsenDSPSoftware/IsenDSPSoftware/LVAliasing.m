function LVAliasing(SR,Freq)
%
% function DemoAliasing(SR,Freq)
%
% This program will generate a
% waveform with a chosen frequency Freq
% and a specified sampling rate SR.
% By repeating the experiment with
% different frequencies for a given
% sampling rate, it will be possible
% to demonstrate the periodic nature of
% the response of the sampling function.
%
%If you choose a sampling rate of 100, try
%various frequencies such as 2, 98, 102, 1002,
%1098, 48, 52, etc.
%
% Sample Calls:
%
% LVAliasing(100,2)
% LVAliasing(100,102)
% LVAliasing(100,98)
% LVAliasing(100,1002)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

t = 0:1/SR:1-1/SR;

Waveform = sin(2*pi*t*Freq);
wavchar = 'Sine';

figure(151)

stem(Waveform,'b.')
TheLim = max(abs(Waveform));

xlabel(['Sample Number'])
ylabel(['Amplitude'])
axis([0 inf -1.1*TheLim 1.1*TheLim])

