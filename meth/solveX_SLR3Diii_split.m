function [x]=solveX_SLR3Diii_b0(x0,k_obj,Ttr,S,A,kGrid,rkGrid,gpu,M)
% 
% Initialise x, z, y;
nX=1;
SH=conj(S);
[kx,ky,kz,coil,shot]=size(x0);
NSte=shot;
NStr=size(A,5);
%yEnd=zeros(kx,ky,kz,1,NSte);
if gpu
    x=gpuArray(x0);A=gpuArray(A);k_obj=gpuArray(k_obj);S=gpuArray(S);SH=gpuArray(SH);
    
M=gpuArray(M);
else
    x=x0;
%     y=par.y;
end


%Hx= Hankel_hb3D(fftdim(x,par.dim(1:2)), par.f,gpu,par.dim); %number of patches * patches * slices
% g=Hx;

%z=Hx*0;

%nz=size(z,3);
% H'H is just voxel-wise scaling by N
% i.e., (H'*H)*x = N.*x
% N  =   adj_hankel_hb3D(ones(size(Hx)), par.f, kx, ky, kz,coil,shot);

% Define shortcut functions for hankel and adjoint
%H_fwd   =   @(x)Hankel_hb3D(fftdim(x,par.dim(1:2)), par.f,gpu,par.dim);
%H_adj   =   @(x)ifftdim(adj_hankel_hb3D(x, par.f, kx, ky, kz,coil,shot,gpu,par.dim),par.dim(1:2));

% Define motion
  
    etDirtr=precomputationsSincRigidTransform(kGrid,[],rkGrid,Ttr,1,0);
    etInvtr=precomputationsSincRigidTransform(kGrid,[],rkGrid,Ttr,0,0);  
%    etDirte=precomputationsSincRigidTransform(kGrid,[],rkGrid,Tte,1,0);
%    etInvte=precomputationsSincRigidTransform(kGrid,[],rkGrid,Tte,0,0);      
if gpu
      	etDirtr{1}=gpuArray(etDirtr{1});
        for m=1:3
            for n=2:3
                etDirtr{n}{m}=gpuArray(etDirtr{n}{m});
            end
        end
	etInvtr{1}=gpuArray(etInvtr{1});
        for m=1:3
            for n=2:3
                etInvtr{n}{m}=gpuArray(etInvtr{n}{m});
            end
        end
       
end
    
    %A'y
    yS=bsxfun(@times,k_obj,A);
    %F'A'y 
    yS=ifftdim(yS,1:3);
    %S'F'A'y
    yS=sum(bsxfun(@times,yS,SH),4);  
    %T'S'F'A'y

        yS=sincRigidTransform_XC(yS,etInvtr,0,gpu);
        for ii=1:NSte
        yEnd(:,:,:,1,ii)=sum(yS(:,:,:,:,(ii-1)*NStr/NSte+1:ii*NStr/NSte),5);
        end
       % yEnd=sincRigidTransform_XC(yEnd,etInvte,0,gpu);
yEnd=bsxfun(@times,yEnd,M);

for i = 1:nX


    
    
    % z subproblem
    % argmin_z { lambda||z||_* + (rho/2)||z - H(x) + y||_2^2
    
    % Solve via singular value soft thresholding
 %   for j=1:nz
 %   [u,s,v] =   svd(Hx(:,:,j) - y(:,:,j), 'econ');
 %   s      =   diag(max(diag(s) - par.lambda/par.rho, 0));
 %   z(:,:,j)=u*s*v';
 %   end
    
%     z_hat   =   par.r*z + (1-par.r)*g;
    
    
    
%     x subproblem
%     argmin_x { (0.5)||Mx-k||_2^2 + (rho/2)||z - H(x) + y||_2^2 }
        % (M'M + (rho)H'H)x = M'k + (rho)H'(z+y)
        % solve directly
    yy   =   reshape((1/2)*yEnd,[],1); 
% pcg
    [k,~]   =   pcg(@(x)cgfun(x), yy, 1E-6, 100,[],[],x(:)); 
    x       =   reshape(k, kx, ky, kz, coil,[]);


   % Hx      =   H_fwd(x);    

    
    % dual update
   % y       =   y + z - Hx;
 
    
    % penalty adjustment
%     if up==1
%    s  = par.rho*H_fwd(x-x0);
%    if norm(z(:)-Hx(:)) > 10*norm(s(:))
%        par.rho = par.rho*par.m;
%        y       = y/par.m;
%    elseif norm(s(:)) > 10*norm(z(:)-Hx(:))
%        par.rho = par.rho/par.m;
%        y       = y*par.m;
%     else
%        up=0;
%    end

%     end
    
  
%    x0 = x;
%     g  = Hx;
    
              
         
end

if gpu
x=gather(x);
end




function q = cgfun(x)
% S*F*Sens*F'
x = reshape(x, kx, ky, kz, coil,[]);

        %Tx
        %xSS=sincRigidTransform_XC(x,etDirte,1,gpu);
        for ii=1:NSte
        xS(:,:,:,:,(ii-1)*NStr/NSte+1:ii*NStr/NSte)=repmat(x(:,:,:,:,ii),[1 1 1 1 NStr/NSte]);
        end
        xS=sincRigidTransform_XC(xS,etDirtr,1,gpu);
        %STx
        xS=bsxfun(@times,xS,S);
        %FSTx
        xS=fftdim(xS,1:3);
        %A'AFSTx
        xS=bsxfun(@times,xS,A);
        %F'A'AFSTx        
        xS=ifftdim(xS,1:3);   
        %S'F'A'AFSTx
        xS=sum(bsxfun(@times,xS,SH),4);
        %T'S'F'A'AFSTx
        xS=sincRigidTransform_XC(xS,etInvtr,0,gpu);
        for ii=1:NSte
        xEnd(:,:,:,1,ii)=sum(xS(:,:,:,:,(ii-1)*NStr/NSte+1:ii*NStr/NSte),5);
        end

        %xEnd=sincRigidTransform_XC(xEnd,etInvte,0,gpu);
xEnd=bsxfun(@times,xEnd,M);
q = (1/2)*xEnd;%+ (par.rho/2)*ifftdim(fftdim(x,par.dim(1:2)).*N,par.dim(1:2));
 
q = reshape(q,[],1); 
end

end


% function cost(i, x, Hx, k_obj, par)
% [~,s,~] = svd(Hx, 0);
% %single coil
% % c1 = 0.5*norm(reshape(x.*par.sample - k_obj, [],1))^2;
% %multiple coils
% c1 = norm(reshape(M_fwd(x,par) - k_obj, [],1))^2;
% c2 = sum(diag(s));
% 
% % c3 = par.L2/2*norm(reshape(x-par.k0,[],1))^2;
% % 
% % kt=permute(reshape(ifftdim(x,1:2),[par.kx*par.ky*par.coilcom,par.shot]),[2,1]);
% % c4=par.TV/2*sum(sum(conj(kt).*(par.L*kt)));
% 
% c = par.w1/2*c1+par.lambda*c2;
% fprintf(1, 'Iter: %d  Cost: %f  c1: %f  c2: %f  Rho: %f\n', i, c, c1,c2,par.rho);
% end

