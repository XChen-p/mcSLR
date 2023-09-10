function [f] = ifft2c(x)

% ifft2c Centric 2D discrete inverse Fourier Transform
%
%    f = ifft2c(x) returns the 2D IFT of matrix X.
%

     f = ifftshift(ifft2(ifftshift(x)));
