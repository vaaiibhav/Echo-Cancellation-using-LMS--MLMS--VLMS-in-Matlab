function LVxZxformFromSamps(b,a,NumzSamps,M,nVals)
% function LVxZxformFromSamps(b,a,NumzSamps,M,nVals)
% Receives a set of coefficients b,a representative of the
% z-transform of a sequence x[n] for which the unit circle is
% included in the ROC. NumzSamps samples of the z-transform of x[n] are
% computed, and then the complete z-transform is computed. From this, nVals 
% samples of x[n] are reconstructed using numerical contour integration.
% The samples of X(z) are also used to reconstruct a periodic version of
% x[n] over nVals samples. Both reconstructed sequences, x[n] and the
% periodic version are plotted for comparison.
% M is number of dense grid z-samples to use for the
% contour integration
% Test call:
% LVxZxformFromSamps([1,1,1,1],[1],6,5000,10)
%
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
dw = 2*pi/(M-1);
RadArg = ([0:dw:2*pi]); 
ZZ = exp(j*RadArg);
dC = ZZ(1,2:length(ZZ)) - ZZ(1,1:length(ZZ)-1); 

k = 0:1:NumzSamps-1; 
Fndx = fix((k*M/NumzSamps)) + 1;
z = ZZ(Fndx);

Nm = 0; Dm = 0;
for nn = 0:-1:-length(b)+1
    Nm = Nm + (b(-nn+1))*(z.^nn);
end
for d = 0:-1:-length(a)+1
    Dm = Dm + (a(-d+1))*(z.^d);
end
zSamps = Nm./Dm;

Nm = 0; Dm = 0; % compute full z-transform
for nn = 0:-1:-length(b)+1
    Nm = Nm + (b(-nn+1))*(ZZ.^nn);
end
for d = 0:-1:-length(a)+1
    Dm = Dm + (a(-d+1))*(ZZ.^d);
end
ZX1 = Nm./Dm;

% regenerate zXform from samples
S = 0;
for k = 0:1:NumzSamps-1
S = S + zSamps(k+1)./((1-(exp(j*2*pi*k/NumzSamps))*(ZZ.^(-1)))+ eps);
end
zXform = (1/NumzSamps)*S.*( 1-ZZ.^(-NumzSamps)+eps ); % recon from samples

zXform(1) = zXform(2) + (zXform(2) - zXform(3)); % avoid singularity at z = 1, i.e., ZZ=0
zXform1 = zXform(1);

figure(4)
subplot(211)
plot(2*[0:1:M-1]/M,abs(ZX1))
xlabel('(a) Frequency, Units of \pi')
ylabel('Magnitude')
axis([0,2,0,1.2*max(abs(ZX1))])

subplot(212)
plot(2*[0:1:M-1]/M,abs(zXform))
xlabel('(b) Frequency, Units of \pi')
ylabel('Magnitude')
axis([0,2,0,1.2*max(abs(zXform))])

avgzXform = 0.5*(zXform(1,2:length(zXform)) + zXform(1,1:length(zXform)-1));
ZZ = ZZ(1:length(zXform));
avgZZ = 0.5*(ZZ(1,2:length(ZZ)) + ZZ(1,1:length(ZZ)-1));
x = zeros(1,nVals);
% do contour integral 
for n = 0:1:nVals-1;
x(n+1) = (1/(2*pi*j))*sum( ((avgZZ.^(n-1)).*avgzXform).*dC); % averaged zXform values 
end

figure(120) % suitable when you know beforehand that x[n] is real

subplot(211)
xvec = 0:1:nVals-1;
stem(xvec,real(x));
ylabel('Real(x[n])')
xlabel('(a) n')

subplot(212) % reconstruction for arbitrary n using direct ifft
lenXtil = length(zSamps);
KK = 0:1:lenXtil-1;
for n = 0:1:nVals-1
    xx(n+1) = real((1/lenXtil)*sum((exp(j*2*pi*n*KK/lenXtil).*zSamps)));
end
stem([0:1:nVals-1],real(xx))
ylabel('Real(x[n])')
xlabel('(b) n')