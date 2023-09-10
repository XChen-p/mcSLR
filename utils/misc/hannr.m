function w = hannr(N)

%   w = hannr(N)
%
%   Symmetric 2D radial Hann window
%   M Chiew
%   Oct 2015

w   =   zeros(N);
for i = 1:N
for j = 1:N
    r       =   sqrt((i-(N/2+1))^2+(j-(N/2+1))^2);
    if r < (N+1)/2
        w(i,j)  = 0.5*(1+cos(2*pi*r/N));   
    end
end
end
