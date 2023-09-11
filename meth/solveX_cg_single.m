function x=solveX_cg_single(x0,k_obj,T,S,A,kGrid,rkGrid,gpu,M)
% 
%k_obj=double(k_obj);
SH=conj(S);
[kx,ky,kz,coil,shot]=size(x0);

if gpu
    x=gpuArray(x0);A=gpuArray(A);k_obj=gpuArray(k_obj);S=gpuArray(S);SH=gpuArray(SH);M=gpuArray(M);
else 
    x=x0;
end





% Define motion
    etDir=precomputationsSincRigidTransform(kGrid,[],rkGrid,T,1,0);
    etInv=precomputationsSincRigidTransform(kGrid,[],rkGrid,T,0,0);   
if gpu
      	etDir{1}=gpuArray(etDir{1});
        for m=1:3
            for n=2:3
                etDir{n}{m}=gpuArray(etDir{n}{m});
            end
        end
	etInv{1}=gpuArray(etInv{1});
        for m=1:3
            for n=2:3
                etInv{n}{m}=gpuArray(etInv{n}{m});
            end
        end
end
    %A'y
    yS=bsxfun(@times,k_obj,A);
    %F'A'y 
    yS=ifftdim(yS,1:3);
%     for m=1:3
%         yS=ifftGPU(yS,m,gpu);
%     end    
    %S'F'A'y
    yS=sum(bsxfun(@times,yS,SH),4);  
    %T'S'F'A'y
    yEnd=sincRigidTransform_XC_single(yS,etInv,0,gpu);   
yEnd=yEnd.*M;
%     yEnd=yS;

    y   =   reshape(yEnd ,[],1); 
% pcg

    [k,~]   =   pcg(@(x)cgfun(x), y, 1E-6, 100,[],[],x(:)); 
    x       =   reshape(k, kx, ky, kz, coil,[]);

if gpu
    x=gather(x);
end
  




function q = cgfun(x)
% S*F*Sens*F'
x = reshape(x, kx, ky, kz, coil,[]);


        %Tx
        xS=sincRigidTransform_XC_single(x,etDir,1,gpu);
% xS=x;
        %STx
        xS=bsxfun(@times,xS,S);
        %FSTx
        xS=fftdim(xS,1:3);
%         for m=1:3
%             xS=fftGPU(xS,m,gpu);
%         end
        %A'AFSTx
        xS=bsxfun(@times,xS,A);
        %F'A'AFSTx        
        xS=ifftdim(xS,1:3);
%         for m=1:3
%             xS=ifftGPU(xS,m,gpu);
%         end   
        %S'F'A'AFSTx
        xS=sum(bsxfun(@times,xS,SH),4);
        %T'S'F'A'AFSTx
        xEnd=sincRigidTransform_XC_single(xS,etInv,0,gpu);
xEnd=xEnd.*M;
% xEnd=xS;


 
q = reshape(xEnd,[],1); 
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

