function [x,par]=solveX_SLR(x0,k_obj,T,S,A,kGrid,rkGrid,nX,gpu,par)
% 
% Initialise x, z, y;
SH=conj(S);
[kx,ky,kz,coil,shot]=size(x0);


Hx= Hankel_pseudo3D(fftdim(x0,1:3), par.f);
g=Hx;

% y = Hx*0;
y=par.y ;
N=par.N;
z=Hx*0;

% H'H is just voxel-wise scaling by N
% i.e., (H'*H)*x = N.*x
%N  =   adj_hankel_pseudo3D(Hankel_pseudo3D(ones(kx,ky,kz,coil,shot), par.f), par.f, kx, ky, kz,coil,shot);

% Define shortcut functions for hankel and adjoint
H_fwd   =   @(x)Hankel_pseudo3D(fftdim(x,1:3), par.f);
H_adj   =   @(x)ifftdim(adj_hankel_pseudo3D(x, par.f, kx, ky, kz,coil, shot),1:3);

% Define motion
    etDir=precomputationsSincRigidTransform(kGrid,[],rkGrid,T,1,0);
    etInv=precomputationsSincRigidTransform(kGrid,[],rkGrid,T,0,0);   
    
    %A'y
    yS=bsxfun(@times,k_obj,A);
    %F'A'y 
    yS=ifftdim(yS,1:3);
    %S'F'A'y
    yS=sum(bsxfun(@times,yS,SH),4);  
    %T'S'F'A'y
    yEnd=sincRigidTransform_XC(yS,etInv,0,gpu);        


for i = 1:nX
    
    % z subproblem
    % argmin_z { lambda||z||_* + (rho/2)||z - H(x) + y||_2^2
    
    % Solve via singular value soft thresholding
    [u,s,v] =   svd(Hx - y, 'econ');
    s      =   diag(max(diag(s) - par.lambda/par.rho, 0));
    z=u*s*v';
    
    z_hat   =   par.r*z + (1-par.r)*g;
    
    
    
%     x subproblem
%     argmin_x { (0.5)||Mx-k||_2^2 + (rho/2)||z - H(x) + y||_2^2 }
        % (M'M + (rho)H'H)x = M'k + (rho)H'(z+y)
        % solve directly
    yy   =   reshape((1/2)*yEnd + (par.rho/2)*H_adj(z_hat+y) ,[],1); 
% pcg
    [k,~]   =   pcg(@(x)cgfun(x), yy, 1E-6, 100); 
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
    
  
    x0 = x;
    g  = Hx;
    
    
%     % print cost
% %     par.verbose=0;
% %     if par.verbose
      
     imm  = sos(reshape(x,[kx,ky*kz,coil*shot]));
%         figure(2)
%       imshow(abs(imm),[])
%       rmse=norm(par.gt(:)-imm(:))/norm(par.gt(:))
%       cost(i, x, Hx, k_obj, par);
% %     end
    
end

par.y=y;

% im  = sos(reshape(x,[par.kx,par.ky,par.coilcom*par.shot])).*par.fac./sqrt(par.state);
% end




function q = cgfun(x)
% S*F*Sens*F'
x = reshape(x, kx, ky, kz, coil,[]);

        %Tx
        xS=sincRigidTransform_XC(x,etDir,1,gpu);
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
        xEnd=sincRigidTransform_XC(xS,etInv,0,gpu);

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

