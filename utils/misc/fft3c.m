function [f] = fft3c(x)

% fft3c Centric 3-dimensional discrete Fourier Transform.
%
%    f = fftnc(x) returns the 3-dimensional DFT of matrix X.
%                 transforms along first 3 dims only.  
%

  n = length(size(x));
  
  if (n <= 3),
    f = fftnc(x);
    
  else,
    f = x;
    for j=1:3,
      f = fftdim(f,j);
    end;
  
  end;
