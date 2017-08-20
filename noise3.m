function sig=noise3(beta, n)
%sig=FractRnd(beta, n)
%To generate 1/f^betta noise
%INPUT
%    beta is the exponent, a number
%    n is the number of synthetic data points
%OUTPUT
%    sig is an n-points array
%
len_beta=length(beta);
phi_n=2*pi*rand(1,n);
%[B,A] = butter(10, 0.2);
%phi_n=filter(B,A,phi_naa);

f=[1:n].^(-beta/2).*exp(i*phi_n(1:n));
sig=real(ifft(f));
sig=(sig-mean(sig))./std(sig);
