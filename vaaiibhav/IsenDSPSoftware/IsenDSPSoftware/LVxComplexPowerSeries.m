function LVxComplexPowerSeries(cn,maxPwr)
% function LVxComplexPowerSeries(cn,maxPwr)
% Raises the complex number cn to the powers
% 0:1:maxPwr and plots the magnitude, the real
% part, the imaginary part, and real v. imaginary parts.
% Test calls:
% LVxComplexPowerSeries(0.69*(1 + j),50)
% LVxComplexPowerSeries(0.99*exp(j*pi/18),40)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
abstotalps = abs(cn.^(0:1:maxPwr));
maxplot = 1.1*max(abstotalps);
minplot = -maxplot;

figure(56)
clf

args = 0:pi/50:2*pi;
cnums = exp(j*args);

subplot(222)
hold on;
plot(real(cnums),imag(cnums),':')
line([0,0],[minplot,maxplot])
line([1.3*minplot,1.3*maxplot],[0,0])

for GrandCtr = 0:1:maxPwr  
ps = cn^GrandCtr;
magPS = abs(ps);

realPS = real(ps);
imagPS = imag(ps);

subplot(221)
hold on
stem(GrandCtr,magPS)
xlabel('(a) Power (i.e., n of exp(jn\theta))')
ylabel('Magnitude')
axis([0,maxPwr,minplot,maxplot])

subplot(222)

hold on
plot(realPS,imagPS,'ro')
xlabel('(b) Real')
ylabel('Imaginary')

axis([1.3*minplot,1.3*maxplot,minplot,maxplot])

subplot(223)
hold on
stem(GrandCtr,realPS)
xlabel('(c) Power')
ylabel('Real')
axis([0,maxPwr,minplot,maxplot])

subplot(224)
hold on
stem(GrandCtr,imagPS)
xlabel('Power')
ylabel('Imaginary')
axis([0,maxPwr,minplot,maxplot])

pause(0.02)
end
