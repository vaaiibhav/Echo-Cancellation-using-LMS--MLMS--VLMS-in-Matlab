function LVxLinearab(a,b,f1,f2,N,Del,LTICoeff)
% function LVxLinearab(a,b,f1,f2,N,Del,LTICoeff)
% Demonstrates the principle of superposition, i.e.,
% LTI(ax1 + by2) = aLTI(x1) + bLTI(x2) where LTI is a
% linear time invariant operator defined by LTICoeff, the
% coefficients of which weight an input sequence and delayed 
% versions thereof, i.e., LTI(x) = c[0]x[n] + c[1]x[n-1] + ..
% where c{n} are the members of the vector LTICoeff.
% a and b are constants, and f1 and f2 are frequencies of cosine waves 
% that are used as x1 and x2 in demonstrating the principle of
% superposition.
% N is the length of the test sequences x1 and x2, and Del is a number of
% samples of delay to impose on x1 and x2 to demonstrate shift invariance.
% Test calls:
% LVxLinearab(2,5,3,5,128,0,[2])
% LVxLinearab(2,-3,13,5,128,0,[2,-1,1,2])
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
figure(78);
clf

t = [0:1:N-1]/N;

x1 = [zeros(1,Del),cos(2*pi*f1*t)];
lenx1 = length(x1);
[y1,nC] = LV_LTIofX(LTICoeff,x1);
leny1 = length(y1);
x2 = [zeros(1,Del),sin(2*pi*f2*t)];
lenx2 = length(x2);
[y2,nC] = LV_LTIofX(LTICoeff,x2);
leny2 = length(y2);

ax1 = a*x1;
bx2 = b*x2;
ay1 = a*y1;
by2 = b*y2;

subplot(421); 
stem(x1,'ko');
xlabel('(a) n')
ylabel('x1')
%axis([ 0 length(x1)  -(1.25*(abs(min(x1))))  (1.25*(abs(max(x1))))  ])
axis([ 0 leny1  -(1.25*(abs(min(y1))))  (1.25*(abs(max(y1))))  ])

subplot(422); 
stem(y1,'ko');
xlabel('(b) n')
ylabel('y1')
axis([ 0 leny1  -(1.25*(abs(min(y1))))  (1.25*(abs(max(y1))))  ])

subplot(423); 
stem(x2,'k*');
xlabel('(c) n')
ylabel('x2')
%axis([ 0 length(x2)  -(1.25*(abs(min(x2))))  (1.25*(abs(max(x2))))  ])
axis([ 0 leny2  -inf  inf ])

subplot(424); 
stem(y2,'ko');
xlabel('(d) n')
ylabel('y2')
axis([ 0 leny2  -(1.25*(abs(min(y2))))  (1.25*(abs(max(y2))))  ])

subplot(425); 
hold on
stem(ax1,'ko');
stem(bx2,'r*');
xlabel('(e) n')
ylabel('ax1, bx2')
axis([ 0 lenx2  -(1.25*(abs(min(min([a*x1,b*x2])))))  (1.25*(abs(max(max([a*x1,b*x2])))))  ])

subplot(426); 
stem(ay1+by2,'ko');
xlabel('(f) n')
ylabel('ay1 + by2')
axis([0 leny2 -(1.25*(abs(min(ay1+by2)))) (1.25*(abs(max(ay1+by2)))) ])

subplot(427); 
stem(ax1+bx2,'ko');
xlabel('(g) n')
ylabel('ax1 + bx2')
axis([ 0 lenx2  -(1.25*(abs(min(ax1+bx2))))  (1.25*(abs(max(ax1+bx2))))  ])

subplot(428); 
[LTIax1Plusbx2,nC] = LV_LTIofX(LTICoeff,(ax1 + bx2));
lenLTIax1Plusbx2 = length(LTIax1Plusbx2);
stem(LTIax1Plusbx2,'ko');
xlabel('(h) n')
ylabel('LTI(ax1 + bx2)')
axis([ 0 lenLTIax1Plusbx2  -(1.25*(abs(min(LTIax1Plusbx2))))  (1.25*(abs(max(LTIax1Plusbx2))))  ])