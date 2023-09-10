function [f] = fft1c(x)

% fft1c Centric one-dimensional discrete Fourier Transform.
%
%    f = fft1c(x) returns the DFT of matrix X.
%

     f = ifftshift(fft(fftshift(x)));
