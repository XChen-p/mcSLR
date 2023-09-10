function M = fftdim_noshift(M,dim)

%
% [m] = fftdim(M,dim)
%
% performs non-centric Fourier Transform along dimension dim.
%

for i = 1:length(dim)
    M   =   fft(M, [], dim(i));
end
