function roc = getROC(tmap, truth)
%
%   roc = getROC(map, truth, [N])
%
%   Mark Chiew
%   Oct 2013
%   Last Updated: Jul 2015
%
%   Produces output receiver-operator characteristic = [FPR, TPR] 
%   Input:
%           map     -   statistical map
%           truth   -   binary truth mask
%           [N]     -   number of bins for threshold parameterisation
%                       (default 100)

N       =   numel(tmap);
thresh  =   sort(tmap(:));
roc_fpr =   zeros(N,1);
roc_tpr =   zeros(N,1);


for i = 1:length(thresh)
    test            =   zeros(size(tmap));
    test(tmap>=thresh(i))    =   1;
    tp      =   nnz(test.*truth);
    fn      =   nnz(test-truth==-1);
    fp      =   nnz(test-truth==1);
    tn      =   nnz(test+truth==0);
    roc_fpr(i)  =   fp/(fp+tn);
    roc_tpr(i)  =   tp/(tp+fn);

end
roc =   [roc_fpr roc_tpr];
