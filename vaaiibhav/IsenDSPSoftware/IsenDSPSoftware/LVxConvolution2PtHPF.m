function LVxConvolution2PtHPF(Freq)
% LVxConvolution2PtHPF(Freq)
% Freq is the frequency of the test sinusoid of length 16 which will be convolved with the 2-point impulse response [1,-1]
% Values of Freq up to 8 will be non-aliased.
% Sample calls: 
% LVConvolution2PtHPF(1)
% LVConvolution2PtHPF(8)
% 
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
 
t = 0:1/16:1-1/16;
yshift(1,1:1:60) = zeros;
ConvolSeqVal(1,1:36) = zeros;
x(1,1:60) = zeros;
y(1,1:60) = zeros;
x(1,41:42) = 1.25*[1 -1];
xabc = 1:1:60;
y(1,25:40)= 0.8*cos(2*pi*t*Freq);

figure(977);

for a = 0:1:24   
   clf
xcv = x(1,25:48);
yshift(1,25:48) = y(1,25-a:48-a);
ConvolSeqVal(1,a+1) = sum(xcv.*yshift(1,25:48));

bb = max(abs(y));
subplot(2,1,1)
hold off
stem(xabc-41,x,'ro');
hold on
stem(xabc-41,yshift,'bo');
xlabel(['(a)  Signal (Amp = 0.8), Time = ',num2str(a-1),'; X-Axis = Samp Time/Seq No.'])
ylabel(['Amp'])
axis([-16 16 -2*bb 2*bb])

subplot(2,1,2)
stem(-1:1:24,ConvolSeqVal(1,1:26),'bo');
aa = 4;
text(4.5,-1.5*aa,['Convolution Sequence Value at Time ',num2str(a-1),' is '...
      ,num2str(ConvolSeqVal(1,a+1))])
xlabel(['(b)  Convolvution Sequence; X-Axis = Sample Time/Sequence Number'])
ylabel(['Amp'])
axis([-1 26 -1.8*aa  1.8*aa])
pause(0.2)
end
