function [psf mip aip]= radial_xf_psf(n_proj,Nr,nt, incr)

if nargin == 3
    incr    =   0;
end

k   =   gen_radial(incr, 2*Nr, n_proj*nt, 1, 360, 1);
k   =   reshape(k,[],nt,2);

psf =   zeros(Nr,Nr,nt);
dims=   [Nr,Nr];

for t = 1:nt
    st  =   nufft_init(reshape(k(:,t,:),[],2),dims,[6,6],dims*2,dims/2);
    w   =   ones(size(k,1),1);
    for ii = 1:20
        tmp =   st.p*(st.p'*w);
        w   =   w./real(tmp);
    end
    psf(:,:,t)  =   nufft_adj(w, st);
end

psf =   fftdim(psf,3);

idx =   setdiff(1:Nr^2,Nr^2/2+Nr/2+1);
tmp =   reshape(abs(psf),[],nt);
mip =   max(tmp,[],1)/max(abs(psf(:)));
aip =   mean(tmp(idx,:),1)/max(abs(psf(:)));
