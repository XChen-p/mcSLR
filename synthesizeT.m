function TGT=synthesizeT(NS,theta,shift,random)

%SYNTHESIZET   Synthesizes a set of rigid transforms affecting the
%   different shots.
%   TGT = SYNTHESIZET(NS,THETA,RANDOMMOTION)
%   synthesizes a set of rigid transforms according to a rotation parameter
%   NS is a vector with each element indicating the number of transforms
%   (corresponding to the number of shots)
%   THETA is a vector with each element specifying the rotation parameter
%   RANDOMMOTION specifies whether motion has to be generated in a random
%   manner by sampling from a uniform distribution in [-theta/2,theta/2]
%   (default) or, provided NS=2 and RANDOMMOTION=0, it should be generated
%   deterministically as two motion states at -theta/2 and theta/2.
%   It returns a cell array of LENGTH(NS) x LENGTH(THETA) transforms
%



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
TGT=bsxfun(@minus,TGT,TGTm);%0-mean due to considerations above ec (23) in the paper
%     end
end
