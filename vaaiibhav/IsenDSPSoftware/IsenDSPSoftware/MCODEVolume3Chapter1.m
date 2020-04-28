% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
% this script is not meant to be called; it contains clean copies of m-code
% from the text. To run any of the example m-code segments below, copy the
% code, paste it into the Command line, and press Return.
return
%Page 10-------------------------------------------------------------------
ML_DragZeros
%Page 11-------------------------------------------------------------------
LVxFreqRespViaZxform([1, 0, -2.0446, 0, 1],1024)
%--------------------------------------------------------------------------
y = abs(fft([1, 0, -2.0446, 0, 1],1024)); plot(y(1,1:512))
%Page 12-------------------------------------------------------------------
a = 0.9*j;b = exp(j*pi/4);Imp = poly([a,-a,(1/a),(1/-a),b,b^-1,-1.0])
%Page 14-------------------------------------------------------------------
LVFrPhRImp101
%Page 19-------------------------------------------------------------------
LVCombFilter(5)
%Page 25-------------------------------------------------------------------
FR = abs(fft([1,zeros(1,24),1],3000)); plot(0:1:180, FR(1,1:181))
%Page 26-------------------------------------------------------------------
[b] = fir1(6,0.5);
[Bc,Ac,Gain] = LVDirToCascade(b,1)
x = chirp([0:1/999:1],0,1,500); y = filter(b,1,x);
[y2] = LVCascadeFormFilter(Bc,Ac,Gain,x);
figure(6); subplot(211); plot(y); subplot(212); plot(y2)
[b,a,k] = LVCas2Dir(Bc,Ac,Gain)
%Page 29-----------------------------------------------------------------
unitImp = [1,zeros(1,6)]; y1 = filter([0.25*[1,0,0,0,-1]],[1],unitImp);
y2 = filter([-1,1],[1,0,1],y1); y3 = filter([3],[1,-1],y1);
y = y2 + y3
