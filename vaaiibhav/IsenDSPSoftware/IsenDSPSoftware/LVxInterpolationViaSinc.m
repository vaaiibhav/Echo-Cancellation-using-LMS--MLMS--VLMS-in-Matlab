function LVxInterpolationViaSinc(N,SampDecRate,valTestSig)
%
% function LVxInterpolationViaSinc(N,SampDecRate,valTestSig)
% N is the master sequence length, from which a densely-sampled test signal
% is generated. The master sequence is decimated by every SampDecRate
% samples to create the test sample sequence from which an interpolated version of
% the underlying bandlimited continuous domain signal will be generated.
% N must be an even integer, and SampDecRate must be an integer.
% valTestSig may be any integer from 1 to 5, with
% 1 yields the waveform cos(2*pi*n*0.1) + 0.7*sin(2*pi*n*0.24);
% 2 yields DC 
% 3 gives a single triangular waveform
% 4 gives two cycles of a triangular waveform
% 5 gives 0.5*cos(2*pi*n*0.125) + 1*sin(2*pi*n*0.22);
% where n = -N/2:1:N/2 if N is even
%
% The master sequence is decimated by the factor SampDecRate to create a sample 
% sequence to be used to reconstruct the(bandlimited) original master sequence using 
% weighted sinc functions.
%
% Sample calls:
%
% LVxInterpolationViaSinc(1000,100,1)
% LVxInterpolationViaSinc(5000,100,1)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
if rem(N,2)==0  % N is even
n = -N/2:1:N/2;
else
error('N must be even')
end

if ~( fix(SampDecRate)-SampDecRate == 0 )
    error('SampDecRate must be an integer')
end

n = n/SampDecRate;
ysinc = sinc(n); 
ysinc = (hamming(length(ysinc))').*ysinc;

%==========================================================================
sincMat([2*N - SampDecRate],fix(N/SampDecRate)) = 0; 
for ctr = 1:1:fix(N/SampDecRate)
sincMat((ctr-1)*SampDecRate + 1:(ctr-1)*SampDecRate + N + 1,ctr) = ysinc';
end
%==========================================================================
szsincMat = size(sincMat);

if valTestSig==1
TestSig = cos(2*pi*n*0.1) + 0.7*sin(2*pi*n*0.24);

elseif valTestSig==2
TestSig = ones(1,length(n));

elseif valTestSig==3
TestSig = [1:1:length(n)]/length(n);
halfLen = fix(length(n)/2);
TestSig = 2*[1:1:halfLen   (halfLen-1):-1:1  ]/length(n);

elseif valTestSig==4
quartLen = fix(length(n)/4);
TestSig = [1:1:quartLen   (quartLen-1):-1:1];
TestSig = [TestSig TestSig];
TestSig = 4*TestSig/length(n);

elseif valTestSig==5
TestSig = 0.5*cos(2*pi*n*0.125) + 1*sin(2*pi*n*0.22);

elseif valTestSig==6
TestSig = 0.33*cos(2*pi*n*0.05) - 0.33*sin(2*pi*n*0.13) + 0.33*sin(2*pi*n*0.29 + pi/4);
end

decTestSig = TestSig(1,1:SampDecRate:length(TestSig));
lendecTestSig = length(decTestSig);
paddedTestSig(1,1:SampDecRate*lendecTestSig) = 0;
paddedTestSig(1,1:SampDecRate:SampDecRate*lendecTestSig) = decTestSig;

figure(222)
clf

decTestSig = decTestSig(1:szsincMat(2));
lendecTestSig = length(decTestSig);

recon = sincMat*decTestSig';

trimrecon = recon(fix(N/2+1):length(recon));

reconsampxvec = [fix(N/2+1):SampDecRate:SampDecRate*lendecTestSig+fix(N/2)];

xvec = 1:SampDecRate:length(TestSig);
xvec = xvec(1,1:lendecTestSig);

hold on
plot(TestSig)
stem(xvec,decTestSig)

plot(trimrecon,'b:')
%grid on
ylabel('Amplitude')
xlabel('Sample')
plotLimHi  = 1.1*max( [ max(decTestSig) max(TestSig)  max(trimrecon) 1 ] );
plotLimLo  = 1.1*min( [ min(decTestSig) min(TestSig)  min(trimrecon)  0 ] );
axis( [0  N  plotLimLo plotLimHi] )


figure(555)
clf

subplot(311)

M = SampDecRate;
hold on
paddedTestSig = zeros(1,M*lendecTestSig);
plot(paddedTestSig)
paddedTestSig(1,1:M:length(paddedTestSig)) = decTestSig;
stem(1:M:length(paddedTestSig), decTestSig);
paddedTestSig(1,1:M:length(paddedTestSig)) = decTestSig;
ylabel('Amplitude')
xlabel('(a) Sample, Zero-Padded Seq')

plotLimHi  = 1.1*max( [ max(paddedTestSig)  1 ] );
plotLimLo  = 1.1*min( [ min(paddedTestSig)   0 ] );
axis( [0  N  plotLimLo plotLimHi] )

ans = conv(ysinc,paddedTestSig);

subplot(312)
plot(ysinc)
ylabel('Amplitude')
xlabel('(b) Sample')
axis([0 length(ysinc)  -0.3 1.1])

subplot(313)
plot(ans)
ylabel('Amplitude')
xlabel('(c) Sample, Convolution of Sinc Fcn & Zero-Padded Sequence')
plotLimHi  = 1.1*max( [ max(ans)  1 ] );
plotLimLo  = 1.1*min( [ min(ans)   0 ] );
axis( [0  length(ans)  plotLimLo plotLimHi] )

figure(333)
clf

hold on
for ctr = 1:1:fix(N/SampDecRate)
plot(sincMat(:,ctr)*decTestSig(ctr)','b:');
end
plot(recon);
ylabel('Amplitude')
xlabel('Sample')

figure(445)
clf

hold on
lim = N/SampDecRate;

szsMat = size(sincMat);

yvec = 0:1:szsMat(1)-1;
for ctr = 1:1:lim
plot(0.5*sincMat(:,ctr) + ctr,yvec,'b.');
end
xlabel(['Sinc Matrix Column'])
ylabel(['Sinc Matrix Row'])
axis([ 0 (1.1*szsMat(2))   0  szsMat(1)])

figure(446)
clf

hold on
for ctr = 1:1:lim
plot(0.5*sincMat(:,ctr)*decTestSig(ctr)'+ ctr,yvec,'b.');
end

if lim < 11
for ctr = 1:1:lim
if ctr-fix(lim/2)<0
    offset = -0.4;
    text(ctr+offset,-80,['x[n ',num2str(ctr-fix(lim/2)),']'])
elseif ctr-fix(lim/2)==0
    offset = -0.25;
    text(ctr+offset,-80,['x[n]'])
else
    offset = -0.4;
    text(ctr+offset,-80,['x[n+',num2str(ctr-fix(lim/2)),']'])
end
end

end
text(-0.1,-80,['\Sigma'])

xlabel(['Sinc Matrix Column'])
ylabel(['Sinc Matrix Row'])
axis([ -1.5 (1.1*szsMat(2))   -99  szsMat(1)])

plot(0.5*recon, yvec);

Lenreconsampxvec = length(reconsampxvec);
lendecTestSig = length(decTestSig);

for ctr = 1:1:lendecTestSig
line([ -1.5   1.05*szsMat(2)  ], [reconsampxvec(ctr)  reconsampxvec(ctr) ]);
end

for rowctr = fix(N/2)+1:SampDecRate:fix(N/2)+1+N-fix(SampDecRate)
matreconaTestSigRow = sincMat(rowctr,:)*decTestSig';
plot(0.5*matreconaTestSigRow,rowctr,'bo')
end

% evaluate error over some interval
diff = ones(1,SampDecRate+1);
for ctr = 1:1:SampDecRate+1
aTestSig = TestSig(N/2  + ctr);
matreconaTestSig = sincMat(N  + ctr,:)*decTestSig';
diff(ctr) = aTestSig - matreconaTestSig;
end
figure(8888)
clf

plot(N:1:N+SampDecRate,diff)
xlabel('Sample of Reconstruction')
ylabel('Error (= Signal - Reconstruction)')
axis([N  N+SampDecRate  -inf  inf])

% stem3 plot of weighted sincs, with sum at highest column number
for ctr = 1:1:lim
wtsincMat(:,ctr) = sincMat(:,ctr)*decTestSig(ctr)';
xMat(:,ctr) = [0:1:szsMat(1)-1]';
yMat(:,ctr) = [ctr*ones(1,szsMat(1))]';
end

wtsincMat(:,lim+1) = sincMat*decTestSig';
xMat(:,lim+1) = [0:1:szsMat(1)-1]';
yMat(:,lim+1) = [(szsMat(2)+1)*ones(1,szsMat(1))]';

figure(222)

return



