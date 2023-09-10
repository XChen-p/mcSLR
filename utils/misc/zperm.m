function z = zperm(data, design, n)

if nargin < 3
    n   =   1000;
end

%   Assume last dimension is time
%   Returns t-statistic map, and dof

dims    =   size(data);
data    =   reshape(double(data),[],dims(end));
%d       =   zeros(length(design),n);
nt      =   dims(2);

%   Step 1 - Filter
%   High-pass filter data using detrend
[B,A]   =   butter(5,0.01,'high');
data    =   filtfilt(B,A,data')';
design  =   filtfilt(B,A,design);

%   Step 2 - Residuals
%   Estimate residuals using OLS
res     =   data - data*pinv(design')*design';

%   Step 3 - Autocorrelation
%   Fit AR(1) model coefficients using Yule-Walker method
ar      =   aryule(res',1);

%   Step 3a - Median filter AR estimates
ar      =   reshape(medfilt2(reshape(ar(:,2),sqrt(dims(1)),sqrt(dims(1))),[3,3]),[],1);

%   Step 4a - Whiten data
%   Apply AR model to data
data    =   data + bsxfun(@times,ar,circshift(data,1,2));
data(:,1)=  data(:,2);

%   Step 4b - Filter design matrix
%   Apply AR model to design
for i = 1:size(design,2)
    d(:,:,i)    =   bsxfun(@plus,design(:,i)',bsxfun(@times,ar,circshift(design(:,i),1)'));
    d(:,1,i)    =   d(:,2,i);
end

%   Loop over voxels

%   Step 5 - Permutations
%   Now that the data are whitened, and exchangability is ensured
%   Randomly permute regressor

%m   =   zeros([size(data,1),1]);
%s   =   zeros([size(data,1),1]);
dd  =   zeros([nt,n]);
z   =   zeros([size(data,1),1]);
    
dd(:,1) =   (1:nt)';
for i = 2:n
    dd(:,i) =   randperm(nt, nt);
end

%{
r   =   var(data,[],2)*nt;
for i = 1:size(design,2)
    d   =   design(:,i);
    c{i}=   data*reshape(d(dd),nt,[]);
    r   =   r - abs(c{i}).^2;
end
%}
%{
for i = 1:size(design,2)
    c{i}=   c{i}./sqrt(r);
end
%}
fprintf(1,'%05d',0);
for x = 1:size(data,1)
if mod(x,16)==0
fprintf(1,'\b\b\b\b\b%05d',x);
end
tmp =   reshape(d(x,dd,:),nt,[],2);
for i = 1:size(design,2)
%{
    tmp =   zeros(nt,n);
    for ii = 1:n
        tmp(:,ii)   =   d{i}(x,dd(:,ii));
    end
    %}
    %tmp =   reshape(d(x,dd,i),nt,[]);
    %   Compute coefficients
    c(i,:)  =   data(x,:)*tmp(:,:,i);
end
r   =   var(data(x,:))*nt;
for i = 1:size(design,2)
    r       =   r - abs(c(i,:)).^2;
end
%r   =   std(r,[],2);
r   =   sqrt(r);
for i = 1:size(design,2)
    %[m s]   =   normfit(c);
    %z(x)    =   (c(1) - m)/s;
    %r       =   std(data(x,:).' - c(x,:).*reshape(d(x,dd), nt, []),[],1);
    [m s]   =   normfit(c(i,:)./r);
    %[m s]   =   normfit(c(i,:));
    z(x,i)  =   (c(i,1)./r(1) - m)/s;
    %z(x,i)  =   (c(i,1) - m)/s;
    %z(x)    =   nnz(c(x,:)>=c(x,1))/n;
end
end
fprintf(1,'\b\b\b\b\b');
