function [I,IB]=sincRigidTransform_XC(I,et,di,gpu)

%SINCRIGIDTRANSFORM   Applies sinc-interpolated 3D rigid transforms
%   [I IB]=SINCRIGIDTRANSFORM(I,ET,DI,GPU) rigid transform of images with 
%   sinc-based interpolation (both forwards and backwards) following M. 
%   Unser, P. Thevenaz, and L. Yaroslavsky, ''Convolution-based 
%   interpolation for fast, high-quality rotation of images'', IEEE
%   Transactions on Image Processing, 4(10):1375–1381, Oct. 1995
%   I is the image to be transformed
%   ET are the diagonal operators as given by
%   precomputationSincRigidTransform
%   DI is a flag to indicate whether to perform direct or inverse transform
%   GPU is a flag that determines whether to use gpu (1) or cpu (0) 
%   computation
%   It returns:
%   I, the transformed image
%   IB, the auxiliary images after different steps
%

tr=[1 3 2;2 1 3];
if di==1             
    %Rotations    
    for m=1:3             
%         I=fftGPU(I,tr(1,m),gpu);
        I=fftdim(I,tr(1,m));
        IB{m}=I;
        I=bsxfun(@times,et{2}{m},I);                
%         I=ifftGPU(I,tr(1,m),gpu);
        I=ifftdim(I,tr(1,m));
%         I=fftGPU(I,tr(2,m),gpu);
        I=fftdim(I,tr(2,m));
        I=bsxfun(@times,I,et{3}{m});
%         I=ifftGPU(I,tr(2,m),gpu);
        I=ifftdim(I,tr(2,m));
%         I=fftGPU(I,tr(1,m),gpu);
        I=fftdim(I,tr(1,m));
        I=bsxfun(@times,I,et{2}{m});
        if m~=3
%             I=ifftGPU(I,tr(1,m),gpu);
            I=ifftdim(I,tr(1,m));
        end
    end    
    %Translation
    for m=1:3     
        if m~=tr(1,3)
%             I=fftGPU(I,m,gpu);
         I=fftdim(I,m);
        end
    end
    IB{4}=I;    
    I=bsxfun(@times,I,et{1});    
%     for m=3:-1:1   
%         I=ifftGPU(I,m,gpu);       
%     end
 I=ifftdim(I,1:3);     
else
    %Back-translation
%     for m=1:3 
%         I=fftGPU(I,m,gpu);    
%     end
    I=fftdim(I,1:3);  
    I=bsxfun(@times,I,et{1});%I've modified this lately    
    for m=3:-1:1  
        if m~=tr(1,3)
%             I=ifftGPU(I,m,gpu);    
            I=ifftdim(I,m);   
        end
    end   
    %Back-rotations
    for m=3:-1:1 
        if m~=3
%             I=fftGPU(I,tr(1,m),gpu);   
         I=fftdim(I,tr(1,m));     
        end
        I=bsxfun(@times,I,et{2}{m});
%         I=ifftGPU(I,tr(1,m),gpu);
        I=ifftdim(I,tr(1,m));
%         I=fftGPU(I,tr(2,m),gpu);
        I=fftdim(I,tr(2,m));
        I=bsxfun(@times,I,et{3}{m});
%         I=ifftGPU(I,tr(2,m),gpu);
        I=ifftdim(I,tr(2,m));  
%         I=fftGPU(I,tr(1,m),gpu);
        I=fftdim(I,tr(1,m));
        I=bsxfun(@times,et{2}{m},I);
%         if tr==1
%         if m==1
%             I=sum(I,5);
%         end
%         end
%         I=ifftGPU(I,tr(1,m),gpu);       
      I=ifftdim(I,tr(1,m));  
    end    
end
