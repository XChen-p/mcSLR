function [t, dof, fstat, dof1, dof2] = tstat(data, design, res)

%   Assume last dimension is time
%   Returns t-statistic map, and dof

dims    =   size(data);
data    =   double(reshape(data,[],dims(end)));
dof     =   size(design,1)-size(design,2);
if nargin == 2
    res =   0;
end

%   Demean data
data    =   detrend(data','constant')';
dof     =   dof-1;

%   Normalise columns of design matrix
for i = 1:size(design,2)
    design(:,i) =   design(:,i)/norm(design(:,i));
end

%   Compute t-stats
t       =   data*pinv(design');
t       =   t./repmat(sqrt(sum(res + abs(data-t*design').^2,2)/dof),1,size(design,2));
t       =   reshape(t,[dims(1:length(dims)-1) size(design,2)]);
t(isnan(t)) =   0;

if nargout > 2
    m   =   data*pinv(design')*design';
    R   =   var(m,[],2)./var(data,[],2);
    
    dof1    =   size(design,2)-1;
    dof2    =   dof;
    fstat   =   reshape((R*dof2)./((1-R)*dof1), [dims(1:length(dims)-1),1]);
end
