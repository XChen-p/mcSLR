function [f] = ifftnc(x)

% ifftnc Centric n-dimensional discrete Fourier Transform.
%
%    f = ifftnc(x) returns the 3-dimensional iDFT of matrix X.
%                 transforms along first 3 dims only.  
%

  f = fftshift(ifftn(ifftshift(x)));
