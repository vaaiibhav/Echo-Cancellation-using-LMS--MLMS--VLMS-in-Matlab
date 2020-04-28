function LVxEvenOddAboutZero(x,n)
% function LVxEvenOddAboutZero(x,n)
% Decomposes the sequence x into even and odd components about zero
% Test call:
% LVEvenOddAboutZero([1 2 3 4],[3 4 5 6])
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
retN = -fliplr(n); % gets indices flipped around n = zero
newNind = min([n,retN]):1:max([n,retN]);
newX = zeros(1,length(newNind));
% determine relative start point of original sequence index relative to new
% minimum index
relStart = min(n)-min(newNind);
relIndXinnewX = relStart+1:relStart+length(x);
newX(1,relIndXinnewX) = x;
xe = (newX + fliplr(newX))/2;
xo = (newX - fliplr(newX))/2;

figure(200)
clf

subplot(311)
plot(newNind,xe,'bo');
for ctr = 1:1:length(newNind)
    line([newNind(ctr)  newNind(ctr)],[0  xe(ctr)])
end
xlabel('(a) Sample Index n')
ylabel('Xe')
axis([ (min(newNind))  (max(newNind))  -(1.5*(abs(min(xe))))  (1.5*(abs(max(xe))))  ])

subplot(312)
plot(newNind,xo,'bo');
for ctr = 1:1:length(newNind)
    line([newNind(ctr)  newNind(ctr)],[0  xo(ctr)])
end
xlabel('(b) Sample Index n')
ylabel('Xo')
axis([ (min(newNind))  (max(newNind))  -(1.5*(abs(min(xo))))  (1.5*(abs(max(xo))))  ])

subplot(313)
plot(newNind,xe+xo,'bo');
for ctr = 1:1:length(newNind)
    line([newNind(ctr)  newNind(ctr)],[0  xo(ctr)+xe(ctr)])
end
xlabel('(c) Sample Index n')
ylabel('Xe + Xo')
axis([ (min(newNind))  (max(newNind))  -(1.5*(abs(min(xe+xo))))  (1.5*(abs(max(xe+xo))))  ])