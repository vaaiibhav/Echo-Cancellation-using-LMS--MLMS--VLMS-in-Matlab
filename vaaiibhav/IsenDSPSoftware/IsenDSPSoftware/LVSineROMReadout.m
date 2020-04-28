function LVSineROMReadout(N,F_ROM,D)
% LVSineROMReadout(32,1,3)
% LVSineROMReadout(32,1,34)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
figure(679)
clf

AddInc = D; % D may be 2,3,4,5,etc.

t = 0:1/N:1-1/N;
y = sin(2*pi*t*F_ROM);

plotlim = max(abs(y));

subplot(211)
hold on
stem(y,'bo');
if F_ROM > 1
xlabel(['(a)  ',num2str(F_ROM),' Cycles of a Sinusoidal Wave over ',num2str(N),' Samples'])
else
xlabel(['(a)  ',num2str(F_ROM),' Cycle of a Sinusoidal Wave over ',num2str(N),' Samples'])   
 end
ylabel(['Amplitude'])
axis([0 N+1 -1.5*plotlim 1.5*plotlim])

DecimatedOutput = zeros(1,length(y));

for Ctr = 0:1:N-1
  
   RawAddress = AddInc*Ctr;
   ActualAddress = mod(RawAddress,N);
   SampleToExtract = y(ActualAddress + 1);
   DecimatedOutput(1,Ctr+1) = SampleToExtract;
  
   subplot(212)
   hold on
    plot(Ctr+1,SampleToExtract,'bo');
    line([Ctr+1  Ctr+1],[0  SampleToExtract])
 
   theMod = mod(AddInc,10);
   
xlabel(['(b)  Output (decimation by a factor of ',num2str(AddInc),' in modulo-',num2str(N),' fashion)'])      
ylabel(['Amplitude'])
axis([0 N+1 -1.5*plotlim 1.5*plotlim])
pause(0.125)
end


 










