function LVxSineROMNonIntDecInterp(N,F_ROM,D)
% LVxSineROMNonIntDecInterp(64,1,2.5)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool

SR = 64;

figure(21110);
clf

SR = N;

AddressIncrement = D;
AddressIncrement= abs(real(AddressIncrement));

if AddressIncrement==0
   AddressIncrement=1
end

t = 0:1/SR:1-1/SR;
xx = 0:1:SR-1;
y = sin(2*pi*t*F_ROM);
plotlim = max(abs(y));

subplot(211)
stem(y,'k.')

if F_ROM > 1
xlabel(['(a)  ',num2str(F_ROM),' Cycles of a Sinusoidal Wave over ',num2str(SR),' Samples'])
else     
xlabel(['(a)  ',num2str(F_ROM),' Cycle of a Sinusoidal Wave over ',num2str(SR),' Samples'])   
 end
ylabel(['Amp'])
axis([0 SR+1 -1.5*plotlim 1.5*plotlim])

DecimatedOutput = zeros(1,1);
aaa = length(DecimatedOutput);   
Ctr = 0;

 while aaa < SR + 1    
   RawAddress = AddressIncrement*Ctr;
   ActualAddress = mod(RawAddress,SR);
   ActualAddressPlus1 = mod(RawAddress+1,SR);
   FirstSampToExtract = y(floor(ActualAddress) + 1);
   SecSampToExtract = y(floor(ActualAddressPlus1) + 1);
   
   SampleToExtract = FirstSampToExtract + (ActualAddress - floor(ActualAddress))*(SecSampToExtract-FirstSampToExtract);
   
   DecimatedOutput(1,Ctr+1) = SampleToExtract;
   aaa = length(DecimatedOutput);
   theSign = sign(SampleToExtract);
   
   
subplot(211)

hold on
stem(ActualAddress+1,SampleToExtract,'kd');
   
subplot(212)
   hold on
   stem(Ctr+1,SampleToExtract,'k.');
   xlabel(['(b) Waveform at (a), Decimated by a Factor of ',num2str(AddressIncrement),', in Modulo-',num2str(N),' fashion'])      
   ylabel(['Amp'])
   axis([0 SR+1 -1.5*plotlim 1.5*plotlim])
   
   pause(0.1)
Ctr = Ctr + 1;   
aaa = length(DecimatedOutput);
end




