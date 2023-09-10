function hdrfind(hdr,str)

names = fieldnames(hdr);
str     =   lower(str);

idx = [];

for i = 1:length(names)
    if isstruct(hdr.(names{i}))
        hdrfind(hdr.(names{i}),str);
    end
    if contains(lower(names{i}),str)
        idx =   cat(1,idx,i);
    end
end

for i = 1:length(idx)
    fprintf(1, '%s :', names{idx(i)});
    tmp =   hdr.(names{idx(i)});
    if isempty(tmp)
        disp('[]');
    else
        disp(tmp);
    end
end