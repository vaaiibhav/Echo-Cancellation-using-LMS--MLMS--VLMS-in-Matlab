function LVxDTFT_MS(x,SampOffset,FreqOffsetExp,M,R,TimeOpt,FreqOpt)
% function LVxDTFT_MS(x,SampOffset,FreqOffsetExp,M,R,TimeOpt,FreqOpt)
%
% Computes and displays the magnitude, phase, real, and imaginary parts of the DTFT
% of the sequence x, evaluated over M samples, then computes the same for a
% modified version of x that has been shifted by SampOffset samples and
% multiplied by a complex exponential FreqOffsetExp.
% 
% Pass R as 1 to evaluate from 0 to pi radians, or
% Pass R as 2 to evaluate from 0 to 2*pi radians
% 
% Pass TimeOpt as 1 to let n (the time indices for x1 and x2) be computed
% as n = -(N-1)/2:1:(N-1)/2 for N odd or
% n = -N/2+1:1:N/2; for even N.
% Pass FreqOpt as 1 for symmetrical frequency computation and display (-R*pi
% to +R*pi, for example) or
% Pass FreqOpt as 2 for frequency comp. and display from 0 to R*pi
% 
% Sample calls:
%
% LVxDTFT_MS([cos(2*pi*25*(0:1:100)/100)],0,1,500,2,1,1) % no sample or freq shift
% LVxDTFT_MS([cos(2*pi*25*(0:1:100)/100)],0,exp(j*2*pi*10*(0:1:100)/100),500,2,1,1)
% LVxDTFT_MS([cos(2*pi*25*(0:1:100)/100)],0,exp(-j*2*pi*10*(0:1:100)/100),500,2,1,1)
% LVxDTFT_MS([cos(2*pi*25*(0:1:100)/100)],0,exp(-j*pi/2),500,2,1,1) % shifts phase
%
% LVxDTFT_MS([exp(j*2*pi*25*(0:1:100)/100)],0,exp(-j*2*pi*10*(0:1:100)/100),500,2,1,1)
% LVxDTFT_MS([cos(2*pi*25*(0:1:100)/100)],0,exp(j*2*pi*12.5*(0:1:100)/100),1000,2,1,1)
% LVxDTFT_MS([1 0 1],2,1,300,2,1,1)    % shift x by 2 samples, no frequency shift
%
% LVxDTFT_MS([1 0 1],0,exp(j*2*pi*1*(0:1:2)/3),100,2,1,1) 
% LVxDTFT_MS([cos(2*pi*25*(0:1:100)/100)],0,exp(j*2*pi*12.5*(0:1:100)/100),200,2,1,1)
%
% LVxDTFT_MS([1 0 1],0,exp(j*2*pi*(0:1:2)/3),300,2,1,1)   % no sample shift, but a
% frequency shift of 2*pi/3 radians
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

figure(43)

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

if ~(TimeOpt==1|TimeOpt==2)
    TimeOpt=1;
end

N = length(x);

if TimeOpt==1
    if rem(N,2)>0  % N is odd
        n = -(N-1)/2:1:(N-1)/2;
    else
        n = -N/2+1:1:N/2; % N is even
    end
elseif TimeOpt==2
    n = 0:1:N-1;
end

W = exp(-j*R*pi/M); 
k = MVec; 
dMat = W.^(n'*k); d = x*dMat;

subplot(421); 
plot(R*MVec/M,abs(d)); 
grid on
xlabel('Norm Freq, x{_1}[n]');
ylabel('Mag')

subplot(423); 
plot(R*MVec/M,angle(d))
grid on
xlabel('Norm Freq, x{_1}[n]'); 
ylabel('Radians')

subplot(425); 
plot(R*MVec/M,real(d));
grid on
xlabel('Norm Freq, x{_1}[n]'); 
ylabel('Real')

subplot(427); 
plot(R*MVec/M,imag(d));
grid on
xlabel('Norm Freq, x{_1}[n]'); 
ylabel('Imag')

x2 = x.*FreqOffsetExp;
W = exp(-j*R*pi/M); k = MVec;
n = n + SampOffset; dMat = W.^(n'*k); d = x2*dMat;

subplot(422); 
plot(R*MVec/M,abs(d));
grid on
xlabel('Norm Freq, x{_2}[n]'); 
ylabel('Mag')

subplot(424); 
plot(R*MVec/M,angle(d))
grid on
xlabel('Norm Freq, x{_2}[n]'); 
ylabel('Radians')

subplot(426); 
plot(R*MVec/M,real(d)); 
grid on
xlabel('Norm Freq, x{_2}[n]'); 
ylabel('Real')

subplot(428); plot(R*MVec/M,imag(d));
grid on
xlabel('Norm Freq, x{_2}[n]'); 
ylabel('Imag')

