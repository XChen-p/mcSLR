function [U2, V2, M] = fast_ica_concat(U, S, V)

%   If you feed this a set of spatial components (Nvox x Ncomponents)
%   Then you get a mixing matrix M as output
%   If your original data was USV' 
%   Then your ICA components are:
%   U_ica = UM
%   V_ica = M'SV' (weighted temporal components)

addpath('~steve/matlab/FastICA_25');

[~, M, ~] = fastica([real(U);imag(U)]', 'approach', 'symm', 'g', 'tanh', 'epsilon', 1e-11, 'maxNumIterations', 3000, 'lastEig', size(U,2));

U2  =   U*M;
V2  =   M'*S*V';
