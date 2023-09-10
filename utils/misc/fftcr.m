function [x] = fftcr(f)

     s = size(f);
     x = zeros(s);

     for j=1:s(1),
	     x(j,:) = fftshift(fft(fftshift(f(j,:))));
     end;
