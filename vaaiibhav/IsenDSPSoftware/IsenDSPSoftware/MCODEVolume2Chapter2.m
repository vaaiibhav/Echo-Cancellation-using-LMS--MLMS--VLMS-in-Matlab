% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
% this script is not meant to be called; it contains clean copies of m-code
% from the text. To run any of the example m-code segments below, copy the
% code, paste it into the Command line, and press Return.
return
%Page 30-------------------------------------------------------------------
radFreq = [0:2*pi/499:2*pi]; z = exp(j*radFreq);
Zxform= 1./(1-0.7*z.^(-1)); plot(radFreq/pi,abs(Zxform))
%Page 33-------------------------------------------------------------------
N=150; z = 1; zpowseq = z.^(-1:-1:-N); x = 0.9.^(1:1:N);
NumAns = sum(zpowseq.*x)
zXformAns = (0.9/z)/(1-0.9*(1/z))
%--------------------------------------------------------------------------
N= 30; b = 0.9; n = -N:1:-1; z = 0.6; sig = -(b.^n); z2n = z.^(-n);
subplot(2,1,1); stem(n,sig); xlabel('n'); ylabel('Signal Amp')
subplot(2,1,2); stem(n,z2n); xlabel('n'); ylabel('z^{-n}')
zXform= -sum( (b.^n).*(z.^(-n)) )
%Page 38-------------------------------------------------------------------
y = filter([0 1],[1,-2,1],[1,zeros(1,20)])
%Page 39-------------------------------------------------------------------
imp = 0.9.^(0:1:100); sig = [1,-0.9]; td = conv(imp,sig)
%--------------------------------------------------------------------------
x = randn(1,8), y = filter([1,1],[1,-1.2,0.81],x);
ans = filter([1,-1.2, 0.81],[1,1],y)
%Page 42-------------------------------------------------------------------
freq = -pi:0.02:pi; z = exp(j*freq);
FR = abs(1 + z.^(-4));
figure(44); subplot(211);
plot(freq/pi,FR)
xlabel('Frequency, Units of pi')
ylabel('Magnitude')
subplot(212);
zplane([1 0 0 0 1],1)
xlabel('Real Part')
ylabel('Imaginary Part')
%Page 43-------------------------------------------------------------------
m = 0:1:3; angs = exp(j*(pi*(2*m +1)/4)); theZeros = 0.5^(0.25)*angs
%--------------------------------------------------------------------------
theZeros = roots([1,0,0,0,0.5])
%Page 44-------------------------------------------------------------------
rts = roots([1,1,1,1,1,1])
%Page 45-------------------------------------------------------------------
theRoots = [j*0.6,-j*0.6];
NetConv = [1];
for ctr = 1:1:length(theRoots)
NetConv = conv(NetConv,[1 -theRoots(ctr)]);
end
theCoeff = NetConv
%--------------------------------------------------------------------------
conv([1, -(0 + 0.9j)],[1, -(0 - 0.9j)])
%Page 46-------------------------------------------------------------------
a = poly([0.9*j,-0.9*j])
%--------------------------------------------------------------------------
zzs = [-1,-1,-1]; pls = [0.9, 0.9*j,-0.9*j]; Denom = poly(pls);
N = length(pls) + 1;M = length(zzs) +1;Num = poly(zzs);
b0 = Num(1); z = 1; NProd = 1; for NCtr = 1:1:M-1;
NProd = NProd*(abs(z - zzs(NCtr))); end
DProd = 1; for DCtr = 1:1:N-1;
DProd = DProd*(abs(z - pls(DCtr))); end
MagResp = abs(b0)*abs(z^(N-M))*NProd/DProd
AltMag = abs(sum(Num.*(z.^(0:-1:-(M-1))))/sum(Denom.* (z.^(0:-1:-(N-1)))))
%Page 47-------------------------------------------------------------------
pmags = abs(roots([1,-1.9, 0.95]))
%--------------------------------------------------------------------------
p = 0.9; ImpResp = filter(1,[1,-p],[1,zeros(1,44)])
StepResp = filter(1,[1,-p],[ones(1,95)])
%Page 50-------------------------------------------------------------------
N=50; ImpResp = filter([1,1],[1,-0.9],[1,zeros(1,50)])
%--------------------------------------------------------------------------
ImpResp = conv([0.9.^(0:1:29)],[1,1])
%--------------------------------------------------------------------------
ImpRespAlt = filter([1,1],[1,-0.9],[1,zeros(1,29)])
%Page 52-------------------------------------------------------------------
[R,p,Ck] = residuez([0.1,0.5,0.1],[1,-1.2,0.81])
%--------------------------------------------------------------------------
[R,p,Ck] = residuez([0.1,0.5,0.1],[1,-1.2,0.81])
n = 0:1:50; x = R(1)*p(1).^n + R(2)*p(2).^n; x(1) = x(1) + Ck;
altx = filter([0.1,0.5,0.1],[1,-1.2,0.81],[1 zeros(1,50)])
diff = x - altx, hold on; stem(x); stem(altx)
%Page 53-------------------------------------------------------------------
LVxNumInvZxform([1],[1,-1.2, 0.81],10000,1,[0:1:50])
%--------------------------------------------------------------------------
LVxNumInvZxform([0,1],[1,-2,1],10000,1.05,[0:1:50])
%Page 55-------------------------------------------------------------------
[R,p,Ck] = residuez([1],[1,0,0.81])
%Page 56-------------------------------------------------------------------
n = 0:1:20; x = 0.5*(0.9*j).^n + 0.5*(-0.9*j).^n
altx = filter(1,[1,0,0.81],[1,zeros(1,20)])
figure; hold on; stem(n,x,'bo'); stem(n,altx,'b*')
%Page 57-------------------------------------------------------------------
w = pi/6; b0 = 1; a1 = - 0.85; Cs = b0*a1/(a1-exp(j*w));
Ce = b0/(1-a1*exp(-j*w)); n = 0:1:50; figure(11);
ys = Cs*a1.^n; yE = Ce*exp(j*w*n); yZ = ys+yE;
subplot(321); stem(real(ys)); subplot(322); stem(imag(ys));
subplot(323); stem(real(yE)); subplot(324); stem(imag(yE));
subplot(325); stem(real(yZ)); subplot(326); stem(imag(yZ));
%Page 58-------------------------------------------------------------------
n = [0:1:49]; y = exp(j*(pi/6)).^n;
ans = filter(1,[1,0.85],y);
figure(8); subplot(211); stem(n,real(ans));
subplot(212); stem(n,imag(ans))
%Page 60-------------------------------------------------------------------
zarg = -pi:2*pi/500:pi; eaz = exp(-j*zarg);
FrqR = (1 + eaz)./(1 - 0.95*eaz);
plot(zarg/pi,abs(FrqR)); xlabel('Frequency, Units of pi')
%Page 62-----------------------------------------------------------
Imp = [1,1,1,1,1,1,1,1]
%--------------------------------------------------------------------------
SR=200; for A = 1:1:SR+1; z = exp(j*(A-1)*2*pi/SR);
B = 0 :-1: -length(Imp)+1; zvec = [z.^B];
zXform(A) = sum(Imp.*zvec); end
plot(abs(zXform))
z = exp(j*(A-1)*pi/SR);
%Page 63-------------------------------------------------------------------
Imp = [1 1 1]; B=0:1:length(Imp)-1; SR = 256;
pzMat = (2*pi*(0:1:SR)/SR)'*(-B);
ZZMat = exp(j*pzMat); zXform= ZZMat*Imp'; plot(abs(zXform))
%--------------------------------------------------------------------------
LVxFreqRespViaZxform([1,1,1,1,1,1,1,1],512)
%Page 64-------------------------------------------------------------------
b = [1, zeros(1,600), 0.7, zeros(1,600), 0.9]
CoeffVec = [1,zeros(1,50),0.7,zeros(1,50),0.9];
SR = 256; zVec = exp(j*pi*(0:1/SR:1));
nzCoef = find(abs(CoeffVec)>0); num = zeros(1,length(zVec));
for Ctr = 1:1:length(nzCoef )
AnsThisCoeff = CoeffVec(nzCoef(Ctr))*(zVec.^(-nzCoef(Ctr)+1));
num = num + AnsThisCoeff; end
figure(888); plot(abs(num))
%Page 65-------------------------------------------------------------------
CoeffVec = [1, -1.27, 0.81];
SR = 256; zVec = exp(j*pi*(0:1/SR:1));
nzCoef = find(abs(CoeffVec)>0); den = zeros(1,length(zVec));
for Ctr = 1:1:length(nzCoef )
AnsThisCoeff = CoeffVec(nzCoef(Ctr))*(zVec.^(-nzCoef(Ctr)+1));
den = den + AnsThisCoeff; end
den = 1./den;
figure(777); plot(abs(den))
%--------------------------------------------------------------------------
LVxPlotZXformMagCoeff([1,1,1,1,1,1,1,1],[1],1,1,512)
%Page 66-------------------------------------------------------------------
LVxPlotZXformMagCoeff([1,1,1,1,1,1,1,1],[1],30,0.95,512)
%Page 69-------------------------------------------------------------------
N = 2000; zarg = 0:2*pi/(N-1):2*pi;
zVec = exp(j*zarg); DenomVec = 1 - 0.9*( zVec.^(-1));
NetzXform= 1./DenomVec;
figure(10001); plot(zarg/pi,abs(NetzXform))
xlabel('Frequency, Units of pi')
%--------------------------------------------------------------------------
LVxPlotZXformMagCoeff([1],[1,-1.8,0.81],[1],[1],256)
%--------------------------------------------------------------------------
LVxPlotZXformMagCoeff([1],[1,0,0.81],[1],[1],256)
%--------------------------------------------------------------------------
help MLPlotZXformMagCoeff
%Page 72-------------------------------------------------------------------
LVxPlotZTransformMag([],[(0.5 + 0.85*j),(0.5 - 0.85*j)],1,1,256)
%Page 73-------------------------------------------------------------------
ansZD = 1./(1 - 0.9*exp(-j*pi/4))
%--------------------------------------------------------------------------
ansTD = sum((0.9.^(0:1:149)).*(exp(-j*pi/4).^(0:1:149)))
%Page 75-------------------------------------------------------------------
n = 0:1:20; x = 0.9.^n;
y = conv([1, -0.9],x)
stem(y)
%--------------------------------------------------------------------------
n = 0:1:20; x = 1.^n;
y = conv([1, -0.9],x)
stem(y)
%--------------------------------------------------------------------------
stem(filter([1,-0.9],[1],[1.^(0:1:50)]))
%Page 82-------------------------------------------------------------------
[b,a] = butter(7,0.4),
x = chirp([0:1/1023:1],0,1,512);
y1 = filter(b,a,x);
[Bc,Ac,Gain] = LVDirToCascade(b,a),
[y2] = LVCascadeFormFilter(Bc,Ac,Gain,x);
figure(10); subplot(211); plot(y1);
subplot(212); plot(y2)
%--------------------------------------------------------------------------
[b,a] = butter(4,0.5)
[Bc,Ac,Gain] = LVDirToCascade(b,a)
[b,a,k] = LVCas2Dir(Bc,Ac,Gain)
%Page 85-------------------------------------------------------------------
[b,a] = butter(5,0.2)
x = chirp([0:1/1023:1],0,1,512);
y1 = filter(b,a,x);
figure(9); subplot(211); plot(y1);
[Bp,Ap,Cp] = LVxDir2Parallel(b,a);
y2 = LVxFilterParallelForm(Bp,Ap,Cp,x);
subplot(212); plot(y2)
%Page 86-------------------------------------------------------------------
[b,a] = butter(3,0.4)
[Bp,Ap,Cp] = LVxDir2Parallel(b,a)
[b,a] = LVxParallel2Dir(Bp,Ap,Cp)
%--------------------------------------------------------------------------
x = chirp([0:1/1023:1],0,1,512);
b = [2,1,0,1.5]; b1 = b(1); b = b/b1; y = filter(b,1,x);
k = tf2latc(b), [g,h] = latcfilt(k,x);
figure(8); subplot(211); plot(y); xlabel('Sample')
subplot(212); plot(g); xlabel('Sample')
b = latc2tf(k); b = b*b1
%Page 87-------------------------------------------------------------------
k = tf2latc([1,0.5,-0.8])
%Page 88-------------------------------------------------------------------
x = chirp([0:1/1023:1],0,1,512);
a = [1,1.3,0.81]; y = filter(1,a,x);
k = tf2latc(1,a),
[g,h] = latcfilt(k,1,x);
figure(8); subplot(211); plot(y);
subplot(212); plot(g)
%Page 89-------------------------------------------------------------------
x = chirp([0:1/1023:1],0,1,512);
a = [1,0.2619,0.2767];
b = [0.3631,0.7263,0.3631];
y = filter(b,a,x);
[k,c] = tf2latc(b,a),
[g,h] = latcfilt(k,c,x);
figure(8); subplot(211); plot(y);
subplot(212); plot(g);
%--------------------------------------------------------------------------


