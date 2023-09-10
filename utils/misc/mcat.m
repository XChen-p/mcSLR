function mcat(input,dim)

if nargin < 2 || isempty(dim)
    dim = ndims(input);
end

M       =   ceil(sqrt(size(input,dim)));
N       =   ceil(size(input,dim)/ceil(sqrt(size(input,dim))));
input   =   squeeze(permute(input,[1:dim-1 dim+1:ndims(input) dim]));
input(:,:,size(input,3)+1:M*N) = 0;

x = size(input,1);
y = size(input,2);

show(reshape(permute(reshape(permute(input,[1,3,2]),x*M,N,[]),[1,3,2]),x*M,y*N));
