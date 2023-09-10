function y = sos(x,dim)

N   =   ndims(x);
if nargin < 2
    dim = N;
end

y   =   sum(abs(x).^2,dim).^0.5;
y   =   permute(y, [setdiff(1:N, dim) dim]);
