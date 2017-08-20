function showFFTPhaseDiffs(varargin)
%showFFTPhaseDiffs Show effects of phase vector in frequency domain

Args = struct('npts',100,'beta',1);

Args = getOptArgs(varargin,Args);

% make sure we have even number of points
npts2 = floor(Args.npts/2);
npts = npts2 * 2;
numUniquePts = npts2 + 1;
endpt = npts2 - 1;

% generate 1/f^beta
f = 1:npts2;
% scale frequencues so that the amplitude of the 1st freq is 128
fexp = f.^(-Args.beta);
% get random phase between -pi and pi
fphase = 2*pi*(rand(1,npts2)-0.5);
fphase2 = 2*pi*(rand(1,endpt)-0.5);

% generate amplitude vector and phase vector
% set DC to get mean of 0 and drop last point of flipped fexp and fphase to 
% make sure mag and phase are npts in length
mag = [0 fexp fliplr(fexp(1:endpt))];
% make sure DC is positive and take negative of fliplr(fphase) to make 
% sure if it is the complex conjugate of the positive frequencies
phase = [0 fphase -fliplr(fphase(1:endpt))];
phase2 = [0 fphase fphase2];

sig = real(ifft( mag.*exp(i*phase) ));
sig2 = real(ifft( mag.*exp(i*phase2) ));

sigf = fft(sig);
sig2f = fft(sig2);
sigfmag = abs(sigf);
sig2fmag = abs(sig2f);
sigfphase = angle(sigf);
sig2fphase = angle(sig2f);

subplot(3,1,1)
plot(sig,'.-')
hold on
plot(sig2,'r.-')
hold off
legend('Conjugate Phase','Random Phase');
title('Signal in Time Domain')

subplot(3,1,2)
plot(mag(1:numUniquePts),'.-')
hold on
plot(sigfmag(1:numUniquePts),'r.-')
plot(sig2fmag(1:numUniquePts),'g.-')
legend('Original','Conjugate Phase','Random Phase');
title('Amplitude in Frequency Domain')
hold off

subplot(3,1,3)
plot(phase(1:numUniquePts),'.-')
hold on
plot(sigfphase(1:numUniquePts),'r.-')
plot(sig2fphase(1:numUniquePts),'g.-')
legend('Original','Conjugate Phase','Random Phase');
title('Phase in Frequency Domain')
hold off
