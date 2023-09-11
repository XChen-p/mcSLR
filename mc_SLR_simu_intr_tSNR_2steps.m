addpath(genpath('.'));
addpath(sprintf('%s/etc/matlab',getenv('FSLDIR')));
pathOu='./data';%Data path
gpu=0;%0->Use CPU / 1->Use GPU
debug=2;%0/1/2 increasing the amount of debug info provided

%Load synthetic data:
% - Ground truth image xGT of size 128x128
% - Sensitivity maps S of size 128x128x1x32
nameX=sprintf('%s/xGT',pathOu);load(nameX);
NS=48;
NStm=48;
NStr=48;% # motion states
NSte=4;
%volume=1;

load('shot1024_phase.mat')
load('S140.mat')
load('epiGT.mat')
xGT=xGT*1E4;
xGTT=bsxfun(@times,xGT,exp(1j*1*phase(:,:,(volume-1)*NS+1:volume*NS)));
Ry=2;Rz=2;
st='st140.nii.gz'
noise=3E-6*1E4;

xGT=permute(xGTT,[1 2 4 5 3]);%128*128*1*1*2

%%%%%%%%%%%%%%%%%%%%
%%%DATA SYNTHESIS%%%
%%%%%%%%%%%%%%%%%%%%

%Transform generation
name=sprintf('TGT302_v%d.mat',volume)
load(name)
TGT=bsxfun(@minus,TGT,TGT(:,:,:,:,1,:));
t4=squeeze(TGT(:,:,:,:,1:8,4));
%tt4 = interp1(t4.',1:0.45:8);
tt4 = interp1(t4.',1:0.146:8);
t1=squeeze(TGT(:,:,:,:,1:8,1));
%tt1 = interp1(t1.',1:0.45:8);
tt1 = interp1(t1.',1:0.146:8);

theta=[0 0 0];% 5 10 20];%Displacement
shift=[0 0 0];
TGT=synthesizeT(NStm,theta,shift,1);%1*1*1*1*shot*6

TGT(:,:,:,:,:,4)=tt4;
TGT(:,:,:,:,:,1)=tt1;

for ii=1:NStm
    TGTT(1,1,1,1,1+(ii-1)*NS/NStm:ii*NS/NStm,:)=repmat(TGT(:,:,:,:,ii,:),[1 1 1 1 NS/NStm 1]);
end
TGT_save=TGT;
TGT=TGTT;

%Grid generation
[kx,ky,kz,~]=size(xGT);
N=[kx,ky,kz];%Image size
[kGrid,kkGrid,rGrid,rkGrid]=generateGrids(N);

%Encoding generation
%[A]=generateEncodingXC(N,NS,EncMeth);
sample = squeeze(CAIPI_Sampling_xc(0,[1,kx,ky,ky/Rz],8,3,Ry,1,1));

A=permute(sample,[1 2 4 5 3]);

%Coil array compression to accelerate
perc=0.99;%Percentage of energy retained
S=coilArrayCompression(S,[],perc,gpu);%128*128*1*11
%Data generation
[y,xS]=synthesizeY(double(xGT),TGT,S,A,kGrid,rkGrid,noise);
sprintf('class(y&S)')
class(xGT)
class(S)
class(y)

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
xid=4;
    mcinter=1;
    mcintra=1;
ls=0;
mean1=0;
mean2=0;
%Solver parameters
mA=1.5;
toler=1E-3;
nExtern=100;%Maximum number of iterations of the joint method
nX=5;%Number of iterations of CG
nT=5;%Number of iterations of Newton's method
winic=1E4;%initial w in ec. (27) of the paper
par.f=6;
par.lambda  = 6E-2;%lambda;%*100;
par.rho     = 1E-1*par.lambda;%0;%1E-15;%;%/10;
% rho adjustment parameter
par.m = 2;%2;
% over-relaxation parameter
par.r = 1;%1.5;


%Spatial mask (which, if available, could be used to constrain solution).
%Disabled here. To enable, the voxels where the solution is constrained to
%be zero should be zero in the mask
M=ones(size(xGT));


%Initialization
x=zeros(kx,ky,kz,1,NSte);
Ttr=zeros(1,1,1,1,NStr,6);
Ttr_save=zeros(1,1,1,1,NStr,6);
Tte=zeros(1,1,1,1,NSte,6);
Tte_d=zeros(1,1,1,1,NSte,6);
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

par.xS=abs(xS(:,:,:,:,1)/maximNormalize);
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

y=zeros((N(1)-par.f+1)*(N(2)-par.f+1), par.f^2*NSte);
NN  =   adj_hankel_pseudo3D(ones(size(y)), par.f, kx, ky, kz,1,NSte);
size(yIn)
size(A)
size(S)
output1 = MCsense(squeeze(sum(yIn,5)),squeeze(sum(A,5)),squeeze(S))*maximNormalize;

for i=1:NSte
    Ate(:,:,:,1,i)=sum(A(:,:,:,:,(i-1)*NStr/NSte+1:i*NStr/NSte),5);
    yte(:,:,:,:,i)=sum(yIn(:,:,:,:,(i-1)*NStr/NSte+1:i*NStr/NSte),5);
end

load(sprintf('Simu302new_mA12_asense_Nmo48_intra48_volume%0.2d.mat',volume))
Ttr(1,1,1,1,:,:)=TtrEst(:,:,end);
%Ttr=TGT;
%Tte(1,1,1,1,:,:)=TtrEst(1:NStr/NSte:end,:,end);
%for i=1:NSte
%    Ttr(1,1,1,1,(i-1)*NStr/NSte+1:i*NStr/NSte,:)=bsxfun(@minus,TtrEst((i-1)*NStr/NSte+1:i*NStr/NSte,:,end),TtrEst((i-1)*NStr/NSte+1,:,end));
%end
clear xEst errX TtrEst
%x=solveX_cg(x,yte,Tte,S,Ate,kGrid,rkGrid,gpu);
sprintf('start')
class(x)
class(S)
class(yIn)
%Ttr_save=bsxfun(@minus,TGT_save,TGT_save(:,:,:,:,1,:));
id=java.util.UUID.randomUUID;
nn=0;
del=0;
nameExp=sprintf('%s/Simu302new_new2_mA12_useasense_moco11_newTte_w0_1d_x%d_Nmo%d_intra%d_inter%d_b%d_volume%0.2d',pathOu,xid,NStm,NStr,NSte,par.lambda*1E3,volume)
for n=1:nExtern
    %Solve for x
    xant=x;%To check convergence
  
    [x,y,par]=solveX_SLRii(x,yIn,Ttr,Tte,S,A,kGrid,rkGrid,nX,y,NN,gpu,par,M);
class(x)
%Solve for T
    if mcinter==1 && n>del
nn=nn+1
            ss1=sprintf('ss1_%s.nii.gz',id);
            save_avw(abs(x(:,:,:,:,1)),ss1,'f',[1 1 1 1])
          for ii=2:NSte            
            ss2=sprintf('ss2_%s.nii.gz',id);            
            save_avw(abs(x(:,:,:,:,ii)),ss2,'f',[1 1 1 1])
            
%            i2rr_new=sprintf('i2rr_new_%s.mat',id);
            i2rr=sprintf('i2rr_%s_%d.mat',id,ii);            
            strr1=strcat('flirt -in ', " ", ss2,' -ref '," ", ss1,' -omat'," ",i2rr,' -2D');
            [status,cmdout] = system(strr1);
            
            Tm=zeros(2,3);
           % if nn==1
           % status = copyfile(i2rr_new, i2rr);
           % else
           %        strr2=strcat('convert_xfm -omat'," ",i2rr,' -concat'," ",i2rr_new," ",i2rr);
           %        [status,cmdout] = system(strr2);
           % end
            
            strr3=strcat('avscale --allparams'," ",i2rr, " ", st,' | head -n 7 | tail -1');
            [~,cmdout] = system(strr3);
            Tm(1,:)=str2num(cmdout(34:end-2));
       	    strr4=strcat('avscale --allparams'," ",i2rr, " ", st,' | head -n 9 | tail -1');            
            [~,cmdout] = system(strr4);
            Tm(2,:)=str2num(cmdout(29:end-2));


            Tte_d(1,1,1,1,ii,4)=Tm(1,3);
            Tte_d(1,1,1,1,ii,1:2)=-Tm(2,1:2);
          end
Tte=Tte+Tte_d;
class(Tte)
if mean1==1   
Tte=Tte-mean(Tte,5);
end
       
end
        
%if xid==3    
% for ii=1:NSte
%        etD=precomputationsSincRigidTransform(kGrid,kkGrid,rkGrid,Tte(:,:,:,:,ii,:),1,0);
%        xt3(:,:,:,:,ii)=sincRigidTransform_XC(x(:,:,:,:,ii),etD,1,gpu);
% end
%end
%if xid==4
% for ii=1:NSte
%       % etD=precomputationsSincRigidTransform(kGrid,kkGrid,rkGrid,Tte(:,:,:,:,ii,:)-mean(Tte,5),1,0);
%       % xt4(:,:,:,:,ii)=sincRigidTransform_XC(x(:,:,:,:,ii),etD,1,gpu);
%etD=precomputationsSincRigidTransform(kGrid,kkGrid,rkGrid,Tte(:,:,:,:,ii,:),1,0);
%xt4(:,:,:,:,ii)=sincRigidTransform_XC(sum(abs(x),5)/NSte.*x(:,:,:,:,ii)./(abs(x(:,:,:,:,ii)+eps)),etD,1,gpu);
%%xt4(:,:,:,:,ii)=sincRigidTransform_XC(sos(abs(x),5)/sqrt(NSte).*x(:,:,:,:,ii)./(abs(x(:,:,:,:,ii)+eps)),etD,1,gpu);
%end
%end

if mcintra==1 && n>del
       
            flagwprev=flagw;%To check convergence
            xx=[];
for ii=1:NSte
if xid==7
%xx=cat(5,xx,repmat(xS(:,:,:,:,1+(ii-1)*NS/NSte)/par.maximNormalize,[1 1 1 1 NStr/NSte]));
xx=cat(5,xx,repmat(sum(abs(x),5)/NSte.*x(:,:,:,:,ii)./(abs(x(:,:,:,:,ii)+eps)),[1 1 1 1 NStr/NSte]));
if n==1
    if ii==NSte
tem=bsxfun(@times,repmat(Tte(:,:,:,:,ii,:)-Tte(:,:,:,:,ii,:),[1 1 1 1 NStr/NSte,1]),permute(1:NStr/NSte,[3 4 5 1 2])/(NStr/NSte));
    else
tem=bsxfun(@times,repmat(Tte(:,:,:,:,ii+1,:)-Tte(:,:,:,:,ii,:),[1 1 1 1 NStr/NSte,1]),permute(1:NStr/NSte,[3 4 5 1 2])/(NStr/NSte));
    end
    Ttr_save(:,:,:,:,(ii-1)*NStr/NSte+1:ii*NStr/NSte,:)= tem;
end

elseif xid==9
xx=cat(5,xx,repmat(sum(abs(x),5)/NSte.*x(:,:,:,:,ii)./(abs(x(:,:,:,:,ii)+eps)),[1 1 1 1 NStr/NSte]));
Ttr_save(:,:,:,:,(ii-1)*NStr/NSte+1:ii*NStr/NSte,:)=Ttr(:,:,:,:,(ii-1)*NStr/NSte+1:ii*NStr/NSte,:)+repmat(Tte(:,:,:,:,ii,:),[1 1 1 1 NStr/NSte 1]);

%if n==1
%    if ii==NSte
%     tem=bsxfun(@times,repmat(Tte(:,:,:,:,ii,:)-Tte(:,:,:,:,ii,:),[1 1 1 1 NStr/NSte,1]),permute(1:NStr/NSte,[3 4 5 1 2])/(NStr/NSte));
%    else
%    tem=bsxfun(@times,repmat(Tte(:,:,:,:,ii+1,:)-Tte(:,:,:,:,ii,:),[1 1 1 1 NStr/NSte,1]),permute(1:NStr/NSte,[3 4 5 1 2])/(NStr/NSte));
%    end
%    Ttr_save(:,:,:,:,(ii-1)*NStr/NSte+1:ii*NStr/NSte,:)= tem+repmat(Tte(:,:,:,:,ii,:),[1 1 1 1 NStr/NSte,1]);
%else
%Ttr_save(:,:,:,:,(ii-1)*NStr/NSte+1:ii*NStr/NSte,:)=Ttr(:,:,:,:,(ii-1)*NStr/NSte+1:ii*NStr/NSte,:)+repmat(Tte(:,:,:,:,ii,:),[1 1 1 1 NStr/NSte,1]);
%end

elseif xid==8
xx=cat(5,xx,repmat(sos(abs(x),5)/sqrt(NSte).*x(:,:,:,:,ii)./(abs(x(:,:,:,:,ii)+eps)),[1 1 1 1 NStr/NSte]));
elseif xid==1
       xx=cat(5,xx,repmat(x(:,:,:,:,ii),[1 1 1 1 NStr/NSte]));%/par.maximNormalize;
elseif xid==2
     xx=cat(5,xx,repmat(x(:,:,:,:,3-ii),[1 1 1 1 NStr/NSte]));%/par.maximNormalize;
elseif xid==3
	etD=precomputationsSincRigidTransform(kGrid,kkGrid,rkGrid,Tte(:,:,:,:,ii,:),1,0);
        xt3=sincRigidTransform_XC(x(:,:,:,:,ii),etD,1,gpu); 
    xx=cat(5,xx,repmat(xt3,[1 1 1 1 NStr/NSte]));%/par.maximNormalize;
%if n==1
%    if ii==NSte
%tem=bsxfun(@times,repmat(Tte(:,:,:,:,ii,:)-Tte(:,:,:,:,ii,:),[1 1 1 1 NStr/NSte,1]),permute(1:NStr/NSte,[3 4 5 1 2])/(NStr/NSte));
%    else
%tem=bsxfun(@times,repmat(Tte(:,:,:,:,ii+1,:)-Tte(:,:,:,:,ii,:),[1 1 1 1 NStr/NSte,1]),permute(1:NStr/NSte,[3 4 5 1 2])/(NStr/NSte));
%    end
%    Ttr_save(:,:,:,:,(ii-1)*NStr/NSte+1:ii*NStr/NSte,:)= tem;
%end


elseif xid==4
etD=precomputationsSincRigidTransform(kGrid,kkGrid,rkGrid,Tte(:,:,:,:,ii,:),1,0);
xt4=sincRigidTransform_XC(sum(abs(x),5)/NSte.*x(:,:,:,:,ii)./(abs(x(:,:,:,:,ii)+eps)),etD,1,gpu);
     xx=cat(5,xx,repmat(xt4,[1 1 1 1 NStr/NSte]));%/par.maximNormalize;
%if n==1
%    if ii==NSte
%tem=bsxfun(@times,repmat(Tte(:,:,:,:,ii,:)-Tte(:,:,:,:,ii,:),[1 1 1 1 NStr/NSte,1]),permute(1:NStr/NSte,[3 4 5 1 2])/(NStr/NSte));
%    else
%tem=bsxfun(@times,repmat(Tte(:,:,:,:,ii+1,:)-Tte(:,:,:,:,ii,:),[1 1 1 1 NStr/NSte,1]),permute(1:NStr/NSte,[3 4 5 1 2])/(NStr/NSte));
%    end
%    Ttr_save(:,:,:,:,(ii-1)*NStr/NSte+1:ii*NStr/NSte,:)= tem;
%end

elseif xid==5
 xx=cat(5,xx,repmat(sum(x,5)/2,[1 1 1 1 NStr/NSte])); 
else 
 xx=cat(5,xx,repmat(sum(abs(x),5)/2,[1 1 1 1 NStr/NSte]));
end
end
Ttr_save=Ttr;
      [Ttr_save,~,w,flagw]=solveT2doublew(xx,yIn,Ttr_save,S,A,kGrid,kkGrid,rGrid,rkGrid,nT,w,flagw,gpu,0,ls);
%        [Ttr_save,~,w,flagw]=solveT2doublerot(xx,yIn,Ttr_save,S,A,kGrid,kkGrid,rGrid,rkGrid,5,w,flagw,gpu,0,ls);
%squeeze(w)     
%        [Ttr_save,~,w2,flagw2]=solveT2doubletra(xx,yIn,Ttr_save,S,A,kGrid,kkGrid,rGrid,rkGrid,5,w2,flagw2,gpu,0,ls);
%squeeze(w2)              

for ii=1:NSte
if mean2==1
   Tmed=mean(Ttr_save(:,:,:,:,(ii-1)*NStr/NSte+1:ii*NStr/NSte,:),5);
else
   Tmed=Ttr_save(:,:,:,:,(ii-1)*NStr/NSte+1,:);
end
Ttr(:,:,:,:,(ii-1)*NStr/NSte+1:ii*NStr/NSte,:)=bsxfun(@minus,Ttr_save(:,:,:,:,(ii-1)*NStr/NSte+1:ii*NStr/NSte,:),Tmed);  
if ii~=1
Tte(:,:,:,:,ii,:)=Tmed+Tte(:,:,:,:,ii,:);
end
end
%Ttr_save=Ttr;



%            for ii=1:NSte
%                    Tmed=-Ttr_save(:,:,:,:,(ii-1)*NStr/NSte+1,:);
%                    Ttr(:,:,:,:,(ii-1)*NStr/NSte+1,:)=0;
%              for jj=2:NStr/NSte
%                Ttr(1,1,1,1,(ii-1)*NStr/NSte+jj,:)=Tsubtract(Tmed,Ttr_save(:,:,:,:,(ii-1)*NStr/NSte+jj,:));%
%              end
%            end
%class(Ttr)
       w(w<1e-4)=2*w(w<1e-4);w(w>1e16)=w(w>1e16)/2;%To avoid numeric instabilities 
%      w2(w2<1e-4)=2*w2(w2<1e-4);w2(w2>1e16)=w2(w2>1e16)/2;   
    end
        %Measure errors
%        xaux=x(:,:,:,:,1)*maximNormalize-xS(:,:,:,:,1);
%        xS1=xS(:,:,:,:,1);
        xaux=sos(x,5)/sqrt(NSte);        
        %             errX(n)=sum(sum(sum((xaux).*conj(xaux))),5)/(N(1)*N(2));
  %          err=norm(xaux(:)-par.xS(:))/norm(par.xS(:))
  %          errX(n)=err;
 %           errF(n)=errorFit_xc(x,Tte,Ttr,S,A,yIn,kGrid,rkGrid);
%errS(n)=ssim(abs(x(:,:,:,:,1)*maximNormalize),abs(xS(:,:,:,:,1)));
            TteEst(:,:,n)=squeeze(Tte);
            TtrEst(:,:,n)=squeeze(Ttr);        
            xEst(:,:,:,n)=squeeze(x)*maximNormalize;
%       wEst(:,n)=squeeze(w);
%       w2Est(:,n)=squeeze(w2);
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
    %xEst=x*maximNormalize; 
     %TteEst=Tte;
     %TtrEst=Ttr;
 
    %delete(s1,s2,i2r,ss1,ss2,ss2_0,i2rr_new,i2rr)

  
  save(nameExp,'output1','xEst','xS','TtrEst','TteEst','TGT','errX','w');
 

