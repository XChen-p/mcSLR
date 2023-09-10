function [f] = ifft1c(x)

% ifft1c Centric one-dimensional discrete inverse Fourier Transform.
%
%    f = ifft1c(x) returns the iDFT of matrix X.
%

     f = fftshift(ifft(ifftshift(x)));
