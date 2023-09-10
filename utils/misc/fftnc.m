function [f] = fftnc(x)

% fftnc Centric n-dimensional discrete Fourier Transform.
%
%    f = fftnc(x) returns the n-dimensional DFT of matrix X.
%

     f = fftshift(fftn(ifftshift(x)));
