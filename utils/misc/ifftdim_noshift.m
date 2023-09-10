function M = ifftdim_noshift(M,dim)

%
% [m] = ifftdim(M,dim)
%
% performs non-centric Fourier Transform along dimension dim.
%

for i = 1:length(dim)
    M   =   ifft(M, [], dim(i));
end
