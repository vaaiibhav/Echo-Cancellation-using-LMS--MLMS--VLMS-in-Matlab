function LVxFreqTest(TestSignalType,N,UserTestSig,dispFreq)
% LVxFreqTest(TestSignalType,N,UserTestSig,dispFreq)
%
% TestSignalType may be passed as 1-7 as follows:
% 1: TestSig = sin(2*pi*t) + 1.25*cos(4*pi*t + pi/3)+ 0.75*cos(12*pi*t + 11*pi/6)
% 2: TestSig = sin(2*pi*t*2 + pi/6) 
% 3: TestSig = sin(2*pi*t*2.71)
% 4: TestSig = 0.25 + cos(2*pi*t*1.31 + pi/2)  +  sin(2*pi*t*2.71 + pi/6)
% 5: TestSig = sin(2*pi*t) + 1.25*cos(2*pi*2*t)+ 0.75*cos(2*pi*5*t)
% 6: TestSig = sin(2*pi*t*2);
% 7 causes UserTestSig to be used as the test signal;  When TestSignalType is 
% passed as other than 7, UserTestSig may be passed as the empty matrix []
% N is the test signal length and must be at least 2. When TestSignalType
% is passed as 7, the value passed for N is overriden by the length of
% UserTestSig. In this case, N may be passed as the empty matrix [] or an 
% arbitrary number if desired.
% Sample Calls: 
% LVxFreqTest(5,32,[],1)
% LVxFreqTest(4,19,[],1)
% LVxFreqTest(7,11,[cos(2*pi*(2.6)*(0:1:10)/11)],1)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

if TestSignalType==7
    Signal = UserTestSig;
    if isempty(Signal)
        error('Option 7 (User test signal) was selected, but no user test signal was supplied')
        return
    end
    N = length(Signal);
    if N<2
       error('User-=supplied test signal must have length of 2 or greater')
        return
    end
else
  N = abs(real(N));  
end

if N<2
    N=2;
end

if rem(N,2)==0
    limN = N/2;
else
    limN = (N-1)/2;
end

t = 0:1/N:1-1/N; % equivalent to n/N with n = 0:1:N-1
CZLReal = zeros(1,limN+1);
CZLImag = zeros(1,limN+1);

if TestSignalType==1
Signal = 1*sin(2*pi*t) + 1.25*cos(2*pi*2*t + 2*pi*(60/360))+ 0.75*cos(2*pi*6*t + 2*pi*(330/360));
elseif TestSignalType==2
Signal = sin(2*pi*t*2 + pi/6);
elseif TestSignalType==3
Signal = sin(2*pi*t*(2.71));
elseif TestSignalType==4
Signal = 0.25 + cos(2*pi*t*(1.31) + pi/2)  +  sin(2*pi*t*(2.71) + pi/6);
elseif TestSignalType==5
Signal = 1*sin(2*pi*t) + 1.25*cos(2*pi*2*t)+ 0.75*cos(2*pi*5*t);
elseif TestSignalType==6
    Signal = sin(2*pi*t*2);
end

for TestFreq = 0:1:limN;
Testcos = cos(2*pi*t*TestFreq);
Testsin = sin(2*pi*t*TestFreq);
CZLReal(1,TestFreq+1) = sum(Signal.*Testcos);
CZLImag(1,TestFreq+1) = -sum(Signal.*Testsin);
if abs(CZLReal(1,TestFreq+1))<10^(-10)
CZLReal(1,TestFreq+1) = 0;
end
if abs(CZLImag(1,TestFreq+1))<10^(-10)
CZLImag(1,TestFreq+1) = 0;
end
%=================================================================
end

figure(650)

subplot(211)
plot(0:1:N-1,Signal,'ko')
hold on
plot(0:1:N-1,Signal,'k')
plot(0:1:N-1,[cos(2*pi*t*dispFreq)],'k')
stem(0:1:N-1,[cos(2*pi*t*dispFreq)],'k*')
hold off
plotStartText = min([0.1*(N-1),0.75]);
axis([0 (N-1),min([-1.5*abs(min(Signal)),-2.5*max(abs(Signal))]),2.5*max(abs(Signal))])

text(plotStartText,2*max(abs(Signal)),['Correlation of Signal (circles) & Test Cosine (stars) = ',num2str(CZLReal(1,dispFreq+1),3)])
xlabel(['(a)  Test Correlator = Cosine, Frequency = ',num2str(dispFreq),'; X-Axis = Sample'])
ylabel(['Amplitude'])

subplot(212)
stem(0:1:N-1,[-sin(2*pi*t*dispFreq)],'k*');

hold on
plot(0:1:N-1,[-sin(2*pi*t*dispFreq)],'k')
plot(0:1:N-1,Signal,'ko')
plot(0:1:N-1,Signal,'k')
text(plotStartText,2*abs(max(Signal)),['Correlation of Signal (circles) & Test Sine (stars) = ',num2str(CZLImag(1,dispFreq+1),3)])
xlabel(['(b)  Test Correlator = Sine, Frequency = ',num2str(dispFreq),'; X-Axis = Sample'])
ylabel(['Amplitude'])
hold off
axis([0,(N-1),min([ -1.5*abs(min(Signal)),-2.5*abs(max(Signal)) ] ),2.5*abs(max(Signal))])

%==========================================================================

figure(651)

SumReal = zeros(1,length(t)); 
SumImag = zeros(1,length(t)); 
% adjust DFT coefficients for bins 0 and N/2 for reconstruction
CZLReal = 2*CZLReal/N;
CZLImag = 2*CZLImag/N;
CZLReal(1,1) = CZLReal(1,1)/2;
if rem(N,2)==0 % N is even, so N/2 exists
   CZLReal(1,(N/2)+1) = CZLReal(1,(N/2)+1)/2; 
end
subplot(311)
%=================Sample-by-Sample, reconstructed waveform=======================
for tCtr = 1:1:length(t)
    thisWave = 0;
for FreqCtr = 0:1:limN  
thisWave = thisWave + CZLReal(1,FreqCtr+1)*cos(2*pi*t(tCtr)*FreqCtr) - CZLImag(1,FreqCtr+1)*sin(2*pi*t(tCtr)*FreqCtr);
end
waveOut(1,tCtr) = thisWave
end
%==================real (cosine) basis vectors weighted by Re(F[k])==========================
currplotmax = 1;
for Ctr = 0:1:limN
   NewAddReal = (CZLReal(1,Ctr+1)*cos(2*pi*t*(Ctr)) );
   thisplotmax = max(abs(NewAddReal));
   currplotmax = max([currplotmax  thisplotmax]);
   hold on
   plot(0:1:N-1,NewAddReal,'k')
   SumReal = SumReal + NewAddReal;
end
xlabel(['(a)  Weighted, Cosine-Correlated Components of Signal; X-Axis = Sample'])
ylabel(['Amplitude'])

plotlim = currplotmax;
axis([0,N-1,-1.5*plotlim,1.5*plotlim])
hold off

subplot(312)
currplotmax = 1;
%=============imag (sine) basis vectors weighted by Im(F[k])====================================
for Ctr1 = 0:1:limN
   NewAddImag = -(CZLImag(1,Ctr1+1)*sin(2*pi*t*(Ctr1)) );
      thisplotmax = max(abs(NewAddImag));
   currplotmax = max([currplotmax  thisplotmax]);
   hold on
   plot(0:1:N-1,NewAddImag,'k')
   SumImag = SumImag + NewAddImag;
end
xlabel(['(b)  Weighted, Sine-Correlated Components of Signal; X-Axis = Sample'])
ylabel(['Amplitude'])

plotlim = currplotmax;

axis([0,N-1,-1.5*plotlim,1.5*plotlim])
hold off

subplot(313)
plot(0:1:N-1,Signal,'ko')
hold on
totalsynthsig = SumReal + SumImag;
plot(0:1:N-1,totalsynthsig,'k*') % synthesized by harmonic vectors
plot(0:1:N-1,waveOut,'k*') % synthesized one sample after another
hold off
xlabel(['(c) Signal (circles); Synth from Coeffs (stars); X-Axis = Sample'])
ylabel(['Amplitude'])
axis([0,inf,-1.2*max(abs(waveOut)),1.2*max(abs(waveOut))])

figure(652)
plotlimreal = max(abs(CZLReal));
plotlimimag = max(abs(CZLImag));

ThePlotLim = max(abs([plotlimreal plotlimimag]));

subplot(311)

xvec = 0:1:N/2;
thedata = CZLReal;
stem(xvec,thedata,'ko');
xlabel(['(a)  Cosine-Correlated Coefficients, X-axis = Frequency'])
ylabel(['Amplitude'])
% axis([0, inf, -1.2*max(abs(thedata)),1.2*max(abs(thedata))])

subplot(312)
xvec = [0:1:N/2];
thedata = CZLImag;
stem(xvec,thedata,'bo');
xlabel(['(b)  Sine-Correlated Coefficients, X-axis = Frequency'])
ylabel(['Amplitude'])
axis([0, inf, min([-1,-1.2*max(abs(thedata))]),max([1,1.2*max(abs(thedata))]) ])

subplot(313)
xvec = [0:1:N/2];
thedata = (CZLImag.^2 + CZLReal.^2).^0.5;
stem(xvec,thedata,'bo');
xlabel(['(c)  Frequency'])
ylabel(['Magnitude'])
axis([0, inf, 0,1.2*max(abs(thedata))])

