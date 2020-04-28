function [b,a] = LVxProny(Imp,NumA,NumB)
% function [b,a] = LVxProny(Imp,NumA,NumB)
% Uses Prony's Method to Model an impulse response as an IIR having z-transform with 
% NumB numerator coefficients b and NumA denominator coefficients a.
% Sample calls:
% [b, a] = LVxProny([(0.9).^( 0:1:50 )],2,2)
% [b, a] = LVxProny([(0.9).^( 0:1:50 )],3,3)
% [b,a] = cheby1(2,0.5,0.5); Imp = filter(b,a,[1,zeros(1,75)]);[b,a]=LVxProny(Imp,3,3)
% [b,a] = cheby1(2,0.5,0.5); Imp = filter(b,a,[1,zeros(1,25)]);[b,a]=LVxProny(Imp,3,3)
% [b,a] = cheby1(2,0.5,0.5); Imp = filter(b,a,[1,zeros(1,10)]);[b,a]=LVxProny(Imp,3,3)
% [b,a] = cheby1(2,0.5,0.5); Imp = filter(b,a,[1,zeros(1,4)]);[b,a]=LVxProny(Imp,3,3)
% [b,a] = cheby1(2,0.5,0.5); Imp = filter(b,a,[1,zeros(1,3)]);[b,a]=LVxProny(Imp,3,3)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
if NumA < 2
    NumA = 2;
end
if NumB < 2
    NumB = 2;
end

Imp = Imp(:)'; lenImp = length(Imp);
HMat = [];
for ctr=1:1:NumA    
  HMat(1:lenImp,ctr) = [zeros(1,ctr-1),Imp(1,1:lenImp-ctr+1)]';
end
Hlp = HMat(NumB+1:lenImp,2:NumA);
RSVec = -HMat(NumB+1:lenImp,1);
Ag0 = Hlp\RSVec;
a = [1,Ag0']';
Hup = HMat(1:NumB,1:NumA);
b = Hup*a;