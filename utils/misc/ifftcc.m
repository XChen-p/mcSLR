function [x] = ifftcc(f)

     s = size(f);
     x = zeros(s);

     for j=1:s(2),
	     x(:,j) = ifftshift(ifft(ifftshift(f(:,j))));
     end;
	     
