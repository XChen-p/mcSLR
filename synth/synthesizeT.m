function TGT=synthesizeT(NS,theta,shift,random)


TGT=zeros([1 1 1 1 NS 6]);%6th dimension for transform parameters: 1-3 translations / 4-6 rotations

for S=1:NS
    if random==0
    TGT(1,1,1,1,S,4:6)=theta*pi/180/(NS-1)*(S-1);
    TGT(1,1,1,1,S,1:3)=shift/(NS-1)*(S-1);
    else
    TGT(1,1,1,1,S,4:6)=theta*pi/180*(rand-0.5);
    TGT(1,1,1,1,S,1:3)=shift*(rand-0.5);
    end
end

TGTm=mean(TGT,5);
TGT=bsxfun(@minus,TGT,TGTm);
end
