function LVxNLSabXSq(a,b,f1,f2,N,Del,NLCoeff)
% function LVxNLSabXSq(a,b,f1,f2,N,Del,NLCoeff)
% Demonstrates that the principle of superposition is not true for a
% nonlinear system defined by NLSCoeff, the
% coefficients of which weight an input sequence and delayed 
% versions thereof raised to the second power, i.e., 
% NLS(x) = c[0]x^2[n] + c[1]x^2[n-1] + c[2]x^2[n-2] + ...
% where c{n} are the members of the vector NLSCoeff.
% a and b are constants, and f1 and f2 are frequencies of cosine 
% and sine waves, respectively, that are used as x1 and x2 in 
% the superposition test i.e., does NLS(ax1 + bx2) = aNLS(x1)+ bNLS(x2)?
% N is the length of the test sequences x1 and x2, and Del is a number of
% samples of delay to impose on x1 and x2 test for shift invariance.
% Test call:
% LVxNLSabXSq(2,-3,13,5,128,0,[2,-1,1,2])
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

figure(88);
clf

t = [0:1:N-1]/N;
x1 = [zeros(1,Del),cos(2*pi*f1*t)];
lenx1 = length(x1);
[y1,nC] = LVxNLSofXSq(NLCoeff,x1);
leny1 = length(y1);
x2 = [zeros(1,Del),sin(2*pi*f2*t)];
lenx2 = length(x2);
[y2,nC] = LVxNLSofXSq(NLCoeff,x2);
leny2 = length(y2);

ax1 = a*x1;
bx2 = b*x2;
ay1 = a*y1;
by2 = b*y2;

subplot(421); 
stem(x1,'ko');
xlabel('(a) n')
ylabel('x1')
axis([ 0 lenx1  -(1.25*(abs(max(y1))))  (1.25*(abs(max(y1))))  ])

subplot(422); 
stem(y1,'ko');
xlabel('(b) n')
ylabel('y1')
axis([ 0 leny1  -(1.25*(abs(max(y1))))  (1.25*(abs(max(y1))))  ])

subplot(423); 
stem(x2,'k*');
xlabel('(c) n')
ylabel('x2')
axis([ 0 lenx2  -(1.25*(abs(max(y2))))  (1.25*(abs(max(y2))))  ])

subplot(424); 
stem(y2,'ko');
xlabel('(d) n')
ylabel('y2')
axis([ 0 leny2  -(1.25*(abs(max(y2))))  (1.25*(abs(max(y2))))  ])

subplot(425); 
hold on
stem(ax1,'ko');
stem(bx2,'r*');
xlabel('(e) n')
ylabel('ax1, bx2')
axis([ 0 lenx2  -(1.25*(abs(max(max([a*x1,b*x2])))))  (1.25*(abs(max(max([a*x1,b*x2])))))  ])

subplot(426); 
stem(ay1+by2,'ko');
xlabel('(f) n')
ylabel('ay1 + by2')
axis([0 leny2 -(1.25*(abs(max(ay1+by2)))) (1.25*(abs(max(ay1+by2)))) ])

subplot(427); 
hold on
stem(ax1+bx2,'ko');
xlabel('(g) n')
ylabel('ax1 + bx2')
axis([ 0 lenx2  -(1.25*(abs(max(ax1+bx2))))  (1.25*(abs(max(ax1+bx2))))  ])

subplot(428); 
[NLSax1Plusbx2,nC] = LVxNLSofXSq(NLCoeff,(ax1 + bx2));
lenNLSax1Plusbx2 = length(NLSax1Plusbx2);
stem(NLSax1Plusbx2,'ko');
xlabel('(h) n')
ylabel('NLS(ax1 + bx2)')
axis([ 0 lenNLSax1Plusbx2  -(1.25*(abs(max(NLSax1Plusbx2))))  (1.25*(abs(max(NLSax1Plusbx2))))  ])