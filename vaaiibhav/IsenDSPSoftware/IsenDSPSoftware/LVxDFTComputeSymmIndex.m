function LVxDFTComputeSymmIndex
%function LVxDFTComputeSymmIndex
%
%This program demonstrates the computation of the Discrete Fourier
%Transform (DFT) using a direct implementation with symmetric indexing of
%k. This program will compute the DFT bin-by-bin;
%the next bin is computed after you press any key.
%Note with the symmetrical indexing scheme, the
%negative frequencies are to the left of zero, and
%the positive frequencies are to the right of zero.
%Contrast this to the asymmetrically-indexed DFT, in
%which the negative frequencies lie to the right of
%Bin N/2. 
%
%The call:
%
% LVxDFTComputeSymmIndex
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
x = -15:1:16;
t = 0:1/32:1-1/32;
A = 1/32;

testvec = sin(2*pi*t)+ 0.3333*sin(2*pi*t*3)+ 0.2*sin(2*pi*t*5)+...
 0.1428*sin(2*pi*t*7)+ 0.1111*sin(2*pi*t*9)+ 0.0909*sin(2*pi*t*11)+...
 0.0769*sin(2*pi*t*13) + 0.0666*sin(2*pi*t*15);

maxplot = 20/32;

figure(23);

DFTrealbin(1,1:32) = 0;
DFTimagbin(1,1:32) = 0;
njk(1,1:32) = 0;

for harmno = -15:1:16
testimag = -sin(2*pi*t*(harmno));
testreal = cos(2*pi*t*(harmno));
DFTrealbin(harmno+16) = A*sum(testvec.*testreal);
DFTimagbin(harmno+16) = A*sum(testvec.*testimag);

subplot(321)
stem(t*32-15,testreal,'bo');
xlabel(['(a)  Test Cosine-Bin ',num2str(harmno),'; X-Axis = Samp No.'])
ylabel(['Amplitude'])
axis([-16 16 -1.3 1.3])

subplot(323)
stem(t*32-15,testimag,'bo');
ylabel(['Amplitude'])
xlabel(['(c)  Test Sine-Bin ',num2str(harmno),'; X-Axis = Samp No.'])
axis([-16 16 -1.3 1.3])

subplot(325)
stem(t*32-15,testvec,'bo');
ylabel(['Amplitude'])
xlabel(['(e) Test Signal; X-Axis = Samp No.'])
axis([-16 16 -1.3 1.3])

subplot(322)
stem(x,DFTrealbin,'bo');
ylabel(['Real(DFT)'])
xlabel(['(b)  Bin'])
axis([-16 16 -maxplot maxplot])

subplot(324)
stem(x,DFTimagbin,'bo');
ylabel(['Imag(DFT)'])
xlabel(['(d)  Bin'])
axis([-16 16 -maxplot maxplot])

subplot(326)
njk = (DFTrealbin.^2 + DFTimagbin.^2).^0.5;
stem(x,njk,'bo');
ylabel(['DFT Mag'])
xlabel(['(f)  Bin'])
axis([-16 16 0 maxplot])

if harmno==16
   return
else
   pause
end

end




