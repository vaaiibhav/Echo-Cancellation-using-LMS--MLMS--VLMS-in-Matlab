function d = LVxDTFT(x,n,M,R,FreqOpt,FigNo)
% function d = LVxDTFT(x,[],NoSamps,2,2,[])
%
% Computes and displays the magnitude, phase, real, and imaginary parts of the DTFT
% of the sequence x having time indices n, evaluated over M samples. To use
% time indices 0:1:N-1, where N is the length of x, pass n as [].
% 
% To not create a display, pass FigNo as [].
% 
% Pass R as 1 to evaluate over pi radians, or
% Pass R as 2 to evaluate over 2*pi radians
% Use FreqOpt = 1 to do a symmetrical frequency evaluation
% or FreqOpt = 2 for an asymmetrical frequency evaluation 
%
% Sample calls:
%
% d = LVxDTFT([fir1(40,[0.2 0.3])],[-20:1:20],500,2,2,88)
% d = LVxDTFT([fir1(40,[0.2 0.3])],[0:1:40],500,2,1,88)
% d = LVxDTFT([cos(2*pi*25*(0:1:99)/100)],[0:1:99],500,2,1,88) 
% d = LVxDTFT([cos(2*pi*25*(0:1:99)/100)],[0:1:99],500,2,1,88)
% d = LVxDTFT([cos(2*pi*25*(-50:1:50)/100)],[-50:1:50],500,2,1,88)
% d = LVxDTFT([cos(2*pi*5*(0:1:20)/20)],[0:1:20],100,2,1,88)
% d = LVxDTFT([cos(2*pi*5*(-10:1:10)/20)],[-10:1:10],100,2,1,88)
% d = LVxDTFT([cos(2*pi*25*(0:1:100)/100)],[-50:1:50],500,2,1,88)
% d = LVxDTFT([exp(j*2*pi*25*(0:1:99)/100)],[0:1:99],500,2,1,88)
% d = LVxDTFT([cos(2*pi*25*(0:1:99)/100)],[0:1:99],1000,2,1,88)
% d = LVxDTFT([1 0 1],[0:1:2],300,2,1,88)   
% d = LVxDTFT([cos(2*pi*25*(0:1:100)/100)],[0:1:100],200,2,1,88)
% d = LVxDTFT([1 0 1],[0:1:2],300,2,1,88)   % no sample shift, but a
% frequency shift of 2*pi/3 radians
% d = LVxDTFT([1 0 1],[],300,2,2,10)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
if ~(FreqOpt==1|FreqOpt==2)
    FreqOpt=1;
end

if FreqOpt==1
    if rem(M,2)>0  % M is odd
        MVec = -(M-1)/2:1:(M-1)/2;
    else
        MVec = -M/2+1:1:M/2; % M is even
    end
elseif FreqOpt==2
    MVec = 0:1:M-1;
end

N = length(x); 

if isempty(n)
    n = 0:1:N-1;
end

W = exp(-j*R*pi/M); 
k = MVec;
dMat = W.^(n'*k); 

d = x*dMat;

if isempty(FigNo)
    return
end

figure(FigNo);
clf

subplot(411); 
plot(R*MVec/M,abs(d)); 
grid on
xlabel(['Normalized Frequency (Multiples of \pi)']); 
ylabel('Mag')

subplot(412); 
plot(R*MVec/M,angle(d))
grid on
xlabel(['Normalized Frequency (Multiples of \pi)']); 
ylabel('Radians')

subplot(413); 
plot(R*MVec/M,real(d));
grid on
xlabel('Normalized Frequency (Multiples of \pi)'); 
ylabel('Real')

subplot(414); 
plot(R*MVec/M,imag(d));
grid on
xlabel('Normalized Frequency (Multiples of \pi)'); 
ylabel('Imag')


