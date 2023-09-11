function [x,y,par]=solveX_SLRiii(x0,k_obj,Ttr,S,A,kGrid,rkGrid,nX,y,N,gpu,par,M)
% 
% Initialise x, z, y;
SH=conj(S);
[kx,ky,kz,coil,shot]=size(x0);
NSte=shot;
NStr=size(A,5);

x=x0;
Hx= Hankel_pseudo3D(fftdim(x,1:3), par.f);
g=Hx;

% y = Hx*0;
% y=par.y ;


z=Hx*0;


% H'H is just voxel-wise scaling by N
% i.e., (H'*H)*x = N.*x
% N  =   adj_hankel_pseudo3D(ones(size(Hx)), par.f, kx, ky, kz,coil,shot);

% Define shortcut functions for hankel and adjoint
H_fwd   =   @(x)Hankel_pseudo3D(fftdim(x,1:3), par.f);
H_adj   =   @(x)ifftdim(adj_hankel_pseudo3D(x, par.f, kx, ky, kz,coil, shot),1:3);

% Define motion
    etDirtr=precomputationsSincRigidTransform(kGrid,[],rkGrid,Ttr,1,0);
    etInvtr=precomputationsSincRigidTransform(kGrid,[],rkGrid,Ttr,0,0);  
    %etDirte=precomputationsSincRigidTransform(kGrid,[],rkGrid,Tte,1,0);
    %etInvte=precomputationsSincRigidTransform(kGrid,[],rkGrid,Tte,0,0);  
    
    %A'y
    yS=bsxfun(@times,k_obj,A);
    %F'A'y 
    yS=ifftdim(yS,1:3);
    %S'F'A'y
    yS=sum(bsxfun(@times,yS,SH),4);  
    %T'S'F'A'y
        yS=sincRigidTransform_XC(yS,etInvtr,0,gpu);
        for ii=1:NSte
        yEnd(:,:,:,:,ii)=sum(yS(:,:,:,:,(ii-1)*NStr/NSte+1:ii*NStr/NSte),5);
        end
yEnd=bsxfun(@times,yEnd,M);
        %yEnd=sincRigidTransform_XC(yEnd,etInvte,0,gpu);
        
for i = 1:nX
    
    % z subproblem
    % argmin_z { lambda||z||_* + (rho/2)||z - H(x) + y||_2^2
    
    % Solve via singular value soft thresholding
    [u,s,v] =   svd(Hx - y, 'econ');
    s      =   diag(max(diag(s) - par.lambda/par.rho, 0));
%     s(s<par.lambda/par.rho)=0;
    z=u*s*v';
    
    z_hat   =   par.r*z + (1-par.r)*g;
    
    
    
%     x subproblem
%     argmin_x { (0.5)||Mx-k||_2^2 + (rho/2)||z - H(x) + y||_2^2 }
        % (M'M + (rho)H'H)x = M'k + (rho)H'(z+y)
        % solve directly
    yy   =   reshape((1/2)*yEnd + (par.rho/2)*H_adj(z_hat+y) ,[],1); 
% pcg
    [k,~]   =   pcg(@(x)cgfun(x), yy, 1E-6, 100,[],[],x(:)); 
    x       =   reshape(k, kx, ky, kz, coil,[]);


    Hx      =   H_fwd(x);    

    
    % dual update
    y       =   y + z_hat - Hx;
 
    
    % penalty adjustment
    s  = par.rho*H_fwd(x-x0);
    if norm(z(:)-Hx(:)) > 10*norm(s(:))
        par.rho = par.rho*par.m;
        y       = y/par.m;
    elseif norm(s(:)) > 10*norm(z(:)-Hx(:))
        par.rho = par.rho/par.m;
        y       = y*par.m;
    end

% up=1;
%         if up==1
%     s  = par.rho*H_fwd(x-x0);
%     if norm(z(:)-Hx(:)) > 10*norm(s(:))
%         par.rho = par.rho*par.m
%         y       = y/m;
%     elseif norm(s(:)) > 10*norm(z(:)-Hx(:))
%         par.rho = par.rho/m
%         y       = y*par.m;
%     else
%        up=0;
%     end


  
    x0 = x;
    g  = Hx;    

%      figure(2)
%     imshow(abs([x(:,:,:,:,1) x(:,:,:,:,2)]),[])

      x1=x(:,:,:,:,1)*par.maximNormalize;
      xS1=par.xS(:,:,:,:,1);
      norm(x1(:)-xS1(:))/norm(xS1(:))
      
end
    par.y=y;
% par.rho=rho;




% im  = sos(reshape(x,[par.kx,par.ky,par.coilcom*par.shot])).*par.fac./sqrt(par.state);
% end




function q = cgfun(x)
% S*F*Sens*F'
x = reshape(x, kx, ky, kz, coil,[]);

        %Tx
%        xSS=sincRigidTransform_XC(x,etDirte,1,gpu);
        for ii=1:NSte
        xS(:,:,:,:,(ii-1)*NStr/NSte+1:ii*NStr/NSte)=repmat(x(:,:,:,:,ii),[1 1 1 1 NStr/NSte]);
        end
        xS=sincRigidTransform_XC(xS,etDirtr,1,gpu);
%     
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
        xEnd(:,:,:,:,ii)=sum(xS(:,:,:,:,(ii-1)*NStr/NSte+1:ii*NStr/NSte),5);
        end
xEnd=bsxfun(@times,xEnd,M);
        %xEnd=sincRigidTransform_XC(xEnd,etInvte,0,gpu);

q = (1/2)*xEnd + (par.rho/2)*ifftdim(fftdim(x,1:3).*N,1:3);
 
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

