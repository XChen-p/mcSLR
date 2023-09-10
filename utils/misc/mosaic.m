function mosaic(input,dim,crange,cmap)

if nargin < 2 || isempty(dim)
    dim = ndims(input);
end
if nargin < 3
    crange = [min(input(:)) max(input(:))];
end
if nargin < 4
    cmap    =   jet;
end

N   =   size(input,dim);
ii  =   repmat({':'},1,ndims(input));

clf;
for i = 1:N
    ii{dim} =   i;
    splt(ceil(sqrt(N)), ceil(N/ceil(sqrt(N))), i, [0 0]);
    show(squeeze(real(input(ii{:}))), crange,'colormap',cmap)
end
