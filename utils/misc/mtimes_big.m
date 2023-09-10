function out = mtimes_big(out, A, B, xx)

if nargin < 4
    n   =   factor(size(out,1));
    xx  =   prod(n(1:round(length(n/2))));
end

for x = reshape(1:size(out,1), [], xx)
    out(x,:)    =   A(x,:)*B;
end
