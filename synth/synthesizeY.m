function [y,xS_gt]=synthesizeY(xGT,TGT,S,A,kGrid,rkGrid,noise_level)

%SYNTHESIZEY   Synthesizes motion corrupted data
%   Y = SYNTHESIZEY(XGT,TGT,S,A,NS,THETA,ENCMETH,XKGRID,KGRID) synthesizes
%   motion corrupted data according to the forward model in Fig. 1 of the
%   paper,
%   XGT is the ground truth image
%   TGT is a cell array of cell array of #TestedShots x #MotionLevels 
%   ground truth transforms
%   S is the coil-array sensitivity map
%   A is a cell array of #TestedShots x #EncodingMethods containing the 
%   sampling scheme
%   KGRID is a 1x3 cell array with the spectral coordinates along each axis
%   RKGRID is a 2x3 cell array where the first dimension of the cell
%   indexes the shearing orientation and the second indexes the 
%   spatio-spectral axes pair used for that particular shearing (given by 
%   [1 3 2;2 1 3]).
%
% for s=1:length(TGT)
%     for v=1:length(TGT{s})  
        %Tx        
        et=precomputationsSincRigidTransform(kGrid,[],rkGrid,TGT,1,0);%1x3cell
        xS=sincRigidTransform_XC(xGT,et,1,0);
        xS_gt=xS;
        %STx
        xS=bsxfun(@times,xS,S);%128*128*1*1*2
        %FSTx
%         for m=1:2
%             xS=fftGPU(xS,m,0);
%         end
        xS=fftdim(xS,1:3);
        %AFSTx
%         for e=1:length(A{s})
%             y=sum(bsxfun(@times,xS,A),5);
y=bsxfun(@times,xS,A);
rng(100);
noise_i = randn(size(y))*noise_level;
noise_j = randn(size(y))*noise_level;
y=noise_i+1j*noise_j+y;
%         end        
%     end
% end
