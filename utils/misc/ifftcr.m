function [x] = ifftcr(f)

     s = size(f);
     x = zeros(s);

     for j=1:s(1),
	     x(j,:) = ifftshift(ifft(ifftshift(f(j,:))));
     end;
	     
