function mask = poisson_disc(dims, r, s, draw)

% 2D poisson disc sampling with minimum distance between points of r
% Mark Chiew (mark.chiew@ndcn.ox.ac.uk)
%
% O(N) algorithm from R. Bridson "Fast Poisson Disk Sampling in Arbitrary
% Dimensions", https://dl.acm.org/citation.cfm?doid=1278780.1278807
%
% Inputs:
%           dims    -   [Nx, Ny] k-space mask dimensions
%           r       -   Minimum spacing between points
%           s       -   Enables variable density sampling 
%                   -   if scalar, is the s.d. of Gaussian scaled by r
%                   -   if array, should have dimensions of dims,
%                       specifying desired density at each location, scaled
%                       between [0,1]
%           draw    -   true/false, whether to draw sampling mask

if nargin < 4
    draw = false;
end
if nargin < 3
    s = r*ones(dims);
end
if nargin < 2
    r = 8;
end

if isscalar(s)
    s = r*exp(((ii-dims(1)/2-1)^2 + (jj-dims(2)/2-1)^2)/(2*s^2));
else
    s = r./s;
end

% step 0 - initialise
mask = false(dims);
list = [];
k = 25;
%{
maskr = false(dims);
for i = 1:2*round(r)+1
    for j = 1:2*round(r)+1
        if (i-round(r))^2 + (j-round(r))^2 < r^2
            maskr(i,j)=true;
        end
    end
end
idx_radius = find(maskr) - sub2ind(dims,round(r),round(r));

maskr = false(dims);
for i = 1:4*round(r)+1
    for j = 1:4*round(r)+1
        if (i-2*round(r))^2 + (j-2*round(r))^2 > r^2 && (i-2*round(r))^2 + (j-2*round(r))^2 < 4*r^2
            maskr(i,j)=true;
        end
    end
end
idx_annulus = find(maskr) - sub2ind(dims,2*round(r),2*round(r));
%}

if draw
    h = imshow(mask);
end

% step 1 - initial sample
%idx = randperm(prod(dims),1);
idx = sub2ind(dims,dims(1)/2+1,dims(2)/2+1);
mask(idx) = true;
list = cat(1,list,idx);

% step 2 - loop
while ~isempty(list)

    idx = randperm(length(list),1);
    [ii,jj] = ind2sub(dims,list(idx));
    
    R = min(s(list(idx)),min(dims)/2);
    
    maskr = false(dims);
    for i = 1:2*round(R)+1
        for j = 1:2*round(R)+1
            if (i-round(R))^2 + (j-round(R))^2 < R^2
                maskr(i,j)=true;
            end
        end
    end
    idx_radius = find(maskr) - sub2ind(dims,round(R),round(R));

    maskr = false(dims);
    for i = 1:4*round(R)+1
        for j = 1:4*round(R)+1
            if (i-2*round(R))^2 + (j-2*round(R))^2 > R^2 && (i-2*round(R))^2 + (j-2*round(R))^2 < 4*R^2
                maskr(i,j)=true;
            end
        end
    end
    idx_annulus = find(maskr) - sub2ind(dims,2*round(R),2*round(R));   

    idx_tmp = min(max(idx_annulus + list(idx),1),prod(dims));
    idx_tmp = idx_tmp(randperm(length(idx_tmp),min(k,length(idx_tmp))));
    match = true;

    for z = 1:length(idx_tmp)    
        %tmp_match = false(dims);
        %tmp_match(min(max(idx_radius + idx_tmp(z),1),prod(dims))) = true;       
        %if ~nnz(tmp_match.*mask)
        %if isempty(intersect(min(max(idx_radius + idx_tmp(z),1),prod(dims)), idx_mask))
        if ~nnz(mask(min(max(idx_radius + idx_tmp(z),1),prod(dims))))
            mask(idx_tmp(z)) = true;
            list = cat(1,list,idx_tmp(z));
            match = false;
            if draw
                pause(1E-6);
                set(h,'CData',mask);
            end
            break;
        end
    end
    if match
        %list = setdiff(list,list(idx));
        list(idx) = [];
    end
    
end