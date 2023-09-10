function out = rms(in)

out =   sqrt(mean(abs(in).^2,ndims(in)));
