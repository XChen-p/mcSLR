function m = fft_win(M,N)
% function m = fft_win(M,N)


% this version only works for restricted cases
s = size(M);
if (s(1)~=s(2)), 
  fprintf('input matrix M must be square\n'); return;
else, 
    P = s(1); 
end;

R = P/N;
if (rem(R,1)~=0), 
  fprintf('matrix size P must be multiple of input size N\n'); return; 
end;


% downsample matrix M by R
v = 1:R:P;
in = [];
for j=0:R:P-R, in = [in j*P+v]; end;

Mdown = zeros(N,N,R*R);
for j=0:R-1,
  for k=0:R-1,
	Mdown(:,:,k*R + j + 1) = reshape( M(in + k*P + j), [N N 1]);
  end;
end;


% take FFTs of submatrices
mdown = zeros(N,N,R*R);
for j=1:R*R,
  mdown(:,:,j) = fft2(ifftshift( Mdown(:,:,j) ));
end;


% calculate complex exponential matrices
v = 0:N-1;
v = -i*2*pi*v/P;
Mexp = zeros(N,N,R);

for k=0:R-1,
  temp = exp(k*v);
  Mexp(:,:,k+1) = reshape( repmat(temp,N,1), [N,N,1] );
end;


% calculate linear combination
m = zeros(N,N);

if 0,
for j=0:R-1,
  for k=0:R-1,
    mtemp = reshape( mdown(:,:,k*R+j+1), [N N] );
    etemp = reshape( Mexp(:,:,j+1), [N N] );
    mtemp = mtemp.*etemp;
    etemp = reshape( Mexp(:,:,k+1), [N N] );
    mtemp = mtemp.*etemp.';
    m = m + mtemp;
  end;
end;

else,
for j=0:R-1,
  for k=0:R-1,
    mtemp = reshape( mdown(:,:,k*R+j+1), [N N] );

    if (j>=R/2) l = R-j+1; else, l = j+1; end;
    etemp = reshape( Mexp(:,:,l+1), [N N] );
    mtemp = mtemp.*etemp;

    if (k>=R/2) l = R-k+1; else, l = k+1; end;
    etemp = reshape( Mexp(:,:,l), [N N] );
    mtemp = mtemp.*etemp.';

    m = m + mtemp;
  end;
end;

end;

m = fftshift(m);
