function [sp_acorr varargout]= spatial_autocorr(data, scale_multiplier, n_timepoints)

%   MChiew
%   Sept 2016
%
%   data is an [Nx, Ny, Nt] array of time-series data
%   scale_multiplier should be >= 1, dictates grid oversampling factor for the autocorrelation
%   n_timepoints is the number of time-points to use to average (chooses them at random)
%
%   sp_acorr is an [Nx, Ny] array containing the 2D spatial autocorrelation 

if nargin < 2
    scale_multiplier=1;
end
if scale_multiplier < 1
    scale_multiplier=1;
end
if nargin < 3
    n_timepoints = 10;
end
if ndims(data) ~= 3
    error('Wrong dimensions of data');
end

dims    =   ceil([size(data,1)*scale_multiplier, size(data,2)*scale_multiplier, size(data,3)]);

%   Compute and rescale mean
m   =   imresize(mean(data,3), scale_multiplier);

%   Compute autocorrelation time-point-by-time-point
acorr   =   zeros(dims(1)*2-1, dims(2)*2-1, dims(3));
for t = randperm(dims(3), n_timepoints)
    acorr(:,:,t)  =   xcorr2(imresize(data(:,:,t),scale_multiplier)-m);
end

%   Final spatial autocorrelation is temporal mean
sp_acorr    =   mean(acorr,3);

%   Spit out x-centre and y-centre line plots if requested
if nargout > 1
    varargout{1}    =   sp_acorr(dims(1),:)/sp_acorr(dims(1),dims(2));
end
if nargout == 3
    varargout{2}    =   sp_acorr(:,dims(2))/sp_acorr(dims(1),dims(2));
end
