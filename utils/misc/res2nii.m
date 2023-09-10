function res2nii(data, name, dims)

data    =   4096*abs(data)/max(abs(data(:)));
save_avw(data, name, 's', dims);
