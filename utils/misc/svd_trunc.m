function out = svd_trunc(input, r)

dims    =   size(input);
[u,s,v] =   lsvd(reshape(input, [], dims(end)), r);
out     =   reshape(u*s*v',dims);
