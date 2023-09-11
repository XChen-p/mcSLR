addpath(genpath('.'));
addpath(sprintf('%s/etc/matlab',getenv('FSLDIR')));
gpu=0;%0->Use CPU / 1->Use GPU
debug=2;%0/1/2 increasing the amount of debug info provided

%Load synthetic data:
NS=48;
NStm=48;
NStr=48;
NSte=1;
%volume=1;

load('shot1024_phase.mat')
load('S140.mat')
load('epiGT.mat')
%xGT=double(xGT);
xGT=xGT*1E4;
xGTT=bsxfun(@times,xGT,exp(1j*1*phase(:,:,(volume-1)*NS+1:volume*NS)));
Ry=2;Rz=2;
st='st140.nii.gz'
noise=3E-6*1E4;
%S=double(S);


xGT=permute(xGTT,[1 2 4 5 3]);%128*128*1*1*2
%size(xGT)


%%%%%%%%%%%%%%%%%%%%
%%%DATA SYNTHESIS%%%
%%%%%%%%%%%%%%%%%%%%

%Transform generation
name=sprintf('TGT302_v%d.mat',volume)
load(name)
%load('TGT24_302.mat')
TGT=bsxfun(@minus,TGT,TGT(:,:,:,:,1,:));
t4=squeeze(TGT(:,:,:,:,:,4));

tt4 = interp1(t4.',1:0.146:8);
%tt4 = interp1(t4.',1:0.315:16);%16
%tt4 = interp1(t4.',1:0.485:24);%24
t1=squeeze(TGT(:,:,:,:,:,1));
tt1 = interp1(t1.',1:0.146:8);

theta=[0 0 0];% 5 10 20];%Displacement
shift=[0 0 0];
TGT=synthesizeT(NStm,theta,shift,1);%1*1*1*1*shot*6

 TGT(:,:,:,:,:,4)=tt4;
 TGT(:,:,:,:,:,1)=tt1;
%TGT=bsxfun(@minus,TGT,mean(TGT(:,:,:,:,1:end,:),5));


%Grid generation
[kx,ky,kz,~]=size(xGT);
N=[kx,ky,kz];%Image size
[kGrid,kkGrid,rGrid,rkGrid]=generateGrids(N);

%Encoding generation
sample = squeeze(CAIPI_Sampling_xc(0,[1,kx,ky,ky/Rz],8,3,Ry,1,1));

A=permute(sample,[1 2 4 5 3]);

%Coil array compression to accelerate
perc=0.99;%Percentage of energy retained
S=coilArrayCompression(S,[],perc,gpu);%128*128*1*11

%Data generation
[y,xS]=synthesizeY(double(xGT),TGT,S,A,kGrid,rkGrid,noise);

%binning
for i=1:NStr
    AA(:,:,:,:,i)=sum(A(:,:,:,:,(i-1)*NS/NStr+1:i*NS/NStr),5);
    yy(:,:,:,:,i)=sum(y(:,:,:,:,(i-1)*NS/NStr+1:i*NS/NStr),5);
end

A=AA;
y=yy;

%%%%%%%%%%%%%%%%%%%%
%%%RECONSTRUCTION%%%
%%%%%%%%%%%%%%%%%%%%
xid=1;
    mcinter=0;
    mcintra=1;
ls=0;
mean1=0;
mean2=0;
%Solver parameters
mA=1.5;
toler=1E-3;
nExtern=200;%Maximum number of iterations of the joint method
nX=5;%Number of iterations of CG
nT=5;%Number of iterations of Newton's method
winic=1E4;%initial w in ec. (27) of the paper
par.f=6;
par.lambda  = 1;%lambda;%*100;
par.rho     = 1E-1*par.lambda;%0;%1E-15;%;%/10;
% rho adjustment parameter
par.m = 2;%2;
% over-relaxation parameter
par.r = 1;%1.5;


%Spatial mask (which, if available, could be used to constrain solution).
%Disabled here. To enable, the voxels where the solution is constrained to
%be zero should be zero in the mask
M=ones(kx,ky,kz,1,NSte);


%Initialization
x=zeros(kx,ky,kz,1,NSte);
Ttr=zeros(1,1,1,1,NStr,6);
Ttr_save=zeros(1,1,1,1,NStr,6);
Tte=zeros(1,1,1,1,NSte,6);
Tr=[];
Tt=[];
NT=size(Ttr);
w=winic*ones(NT(1:6));
w2=winic*ones(NT(1:6));
flagw=zeros(NT(1:6));
flagw2=zeros(NT(1:6));
flagwprev=zeros(NT(1:6));

%Precomputations
yX=ifftdim(y,1:3);
maximNormalize=max(max(max(max(abs(yX)))));
yIn=y/maximNormalize;

par.xS=abs(xGT(:,:,:,:,1)/maximNormalize);
size(par.xS)
par.maximNormalize=maximNormalize;

%xini=zeros(kx,ky,kz,1,NStr);
%xini=repmat(sum(ifftdim(sum(yIn,5),1:2).*conj(S),4),[1 1 1 1 NStr]);
%par.y=zeros((N(1)-par.f+1)*(N(2)-par.f+1), par.f^2*NStr);
%par.N  =   adj_hankel_pseudo3D(ones(size(par.y)), par.f, kx, ky, kz,1,NStr);
%  [xini,~]=solveX_SLR(xini,yIn,Ttr,S,A,kGrid,rkGrid,nX,gpu,par);
% Tini=zeros(1,1,1,1,NStr,6);
%  for ii=1:NSte
% Tini(:,:,:,:,(ii-1)*NStr/NSte+2:ii*NStr/NSte,:)=MCFLIRT(xini(:,:,:,:,(ii-1)*NStr/NSte+1),xini(:,:,:,:,(ii-1)*NStr/NSte+2:ii*NStr/NSte));
% end
%Ttr=Tini;

%y=zeros((N(1)-par.f+1)*(N(2)-par.f+1), par.f^2*NSte);
%NN  =   adj_hankel_pseudo3D(ones(size(y)), par.f, kx, ky, kz,1,NSte);

%Ttr_save=bsxfun(@minus,TGT_save,TGT_save(:,:,:,:,1,:));
id=java.util.UUID.randomUUID;
nameExp=sprintf('%s/Simu302new_changex_asense_1d_Nmo%d_intra%d_volume%0.2d',pathOu,NStm,NStr,volume)
%nameExp=sprintf('%s/Simu000new_mA12_asense_1d_Nmo%d_intra%d_volume%0.2d',pathOu,NStm,NStr,volume)
for n=1:nExtern
    %Solve for x
    xant=x;%To check convergence
    x=solveX_cg_single(x,yIn,Ttr,S,A,kGrid,rkGrid,gpu,M);  
%Solve for T
        

if mcintra==1
       
            flagwprev=flagw;%To check convergence
       xx=repmat(x,[1 1 1 1 NStr/NSte]));%/par.maximNormalize;
      [Ttr,~,w,flagw]=solveT2doublew(xx,yIn,Ttr,S,A,kGrid,kkGrid,rGrid,rkGrid,nT,w,flagw,gpu,0,ls);
       % [Ttr,~,w,flagw]=solveT2doublerot(xx,yIn,Ttr,S,A,kGrid,kkGrid,rGrid,rkGrid,5,w,flagw,gpu,0,ls);
       % [Ttr,~,w2,flagw2]=solveT2doubletra(xx,yIn,Ttr,S,A,kGrid,kkGrid,rGrid,rkGrid,5,w2,flagw2,gpu,0,ls);
               

for ii=1:NSte
if mean2==1
   Tmed=mean(Ttr(:,:,:,:,(ii-1)*NStr/NSte+1:ii*NStr/NSte,:),5);
else
   Tmed=Ttr(:,:,:,:,(ii-1)*NStr/NSte+1,:);
end
Ttr(:,:,:,:,(ii-1)*NStr/NSte+1:ii*NStr/NSte,:)=bsxfun(@minus,Ttr(:,:,:,:,(ii-1)*NStr/NSte+1:ii*NStr/NSte,:),Tmed);  

 etD=precomputationsSincRigidTransform(kGrid,kkGrid,rkGrid,Tmed,1,0);
 x=sincRigidTransform_XC(x,etD,1,gpu);
end
       w(w<1e-4)=2*w(w<1e-4);w(w>1e16)=w(w>1e16)/2;%To avoid numeric instabilities 
       end
        %Measure errors
        xaux=abs(x);
        
        %             errX(n)=sum(sum(sum((xaux).*conj(xaux))),5)/(N(1)*N(2));
%            err=norm(xaux(:)-par.xS(:))/norm(par.xS(:))
%            errX(n)=err;
%            errF(n)=errorFit_xc(x,Tte,Ttr,S,A,yIn,kGrid,rkGrid);
            TtrEst(:,:,n)=squeeze(Ttr);        
            %wEst(:,:,n)=squeeze(w);
       % xant=x-xant;
       % xant=real(xant.*conj(xant));
       % xant=max(max(max(xant)));
        xant=norm(x(:)-xant(:))/norm(xant(:));
        if debug==2;fprintf('Iteration %04d - Error %0.2g \n',n,xant);end
        if xant<toler && sum(flagwprev(:)~=1)>0
            if debug>0;fprintf('Convergence reached at iteration %04d \n',n);end
            break
        elseif n==nExtern
            if debug>0;fprintf('Maximum number of iterations reached without convergence\n');end
        end
    end
    xEst=x*maximNormalize;
 
   
 
    %delete(s1,s2,i2r,ss1,ss2,ss2_0,i2rr_new,i2rr)

  
  save(nameExp,'xEst','xGT','TtrEst','TGT','errX');
 

