function [f] = fft2c(x)

% fft2c Centric two-dimensional discrete Fourier Transform.
%
%    f = fft2c(x) returns the 2DFT of matrix X.
%

     f = fftshift(fft2(fftshift(x)));
