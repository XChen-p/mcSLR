function [f] = ifft3c(x)

% ifft3c Centric 3-dimensional discrete Fourier Transform.
%
%    f = ifftnc(x) returns the 3-dimensional iDFT of matrix X.
%                 transforms along first 3 dims only.  
%

  n = length(size(x));
  
  if (n <= 3),
    f = ifftnc(x);
    
  else,
    f = x;
    for j=1:3,
      f = ifftdim(f,j);
    end;
  
  end;
