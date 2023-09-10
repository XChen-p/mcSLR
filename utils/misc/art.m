function x = art(A, b, niters, lambda)

%   Algebraic Reconstruction Technique
%   Solves for x, in the linear problem b = Ax
%
%   M Chiew
%   Nov 2015

[m,n]   =   size(A);

if nargin == 2
    niters  =   m*n;
    lambda  =   @(i)1;
end

x   =   zeros(n,1);

for k = 1:niters
    i   =   mod(k,m)+1;
    x   =   x + lambda(k)*((b(i)-A(i,:)*x)/(A(i,:)*A(i,:)'))*A(i,:)';
end
