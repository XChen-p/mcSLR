function [G,GB,GC]=sincRigidTransformGradient_xc(xB,et,etg)

%SINCRIGIDTRANSFORMGRADIENT   Computes the gradient of sinc-interpolated 
%3D rigid transforms
%   [G,GB,GC]=SINCRIGIDTRANSFORMGRADIENT(xB,ET,ETG) obtains the 
%   gradient of the transform of the images
%   XB contains images before the first, second and third rotations and 
%   before the translation
%   ET are the diagonal operators of the transform as given by
%   precomputationSincRigidTransform
%   ETG are the diagonal operators of the transform gradient as given by
%   precomputationSincRigidTransform
%   dim is a flag that determines whether to use dim (1) or cpu (0) 
%   computation
%   It returns:
%   G, the gradient of the transformed image
%   GB, the gradient of the first, second and third rotations before 
%   applying the translation
%   GC, the gradient of the first, first and second rotations before
%   applying the second, third and third rotations respectively
%

GB=[];
GC=[];

%Translation parameters
for m=1:3    
    x{1}=bsxfun(@times,xB{4},etg{1}{m});
%     for n=1:3       
%         x{1}=ifftdim(x{1},n);
%     end    
    x{1}=ifftdim(x{1},1:3);
    G{m}=x{1};
end

%First rotation
x{1}=bsxfun(@times,xB{1},et{2}{1}); %tanX
x{2}=bsxfun(@times,xB{1},etg{2}{1});%tan'x
for m=1:2
    x{m}=ifftdim(x{m},1);%F1'tanX, F1'tan'X
    x{m}=fftdim(x{m},2);%x{1}=F2F1'tanX,x{2}=F2F1'tan'X
end
x{2}=bsxfun(@times,etg{3}{1},x{1})+bsxfun(@times,et{3}{1},x{2});%sin'F2F1HtanX + sinF2F1Htan'X
x{1}=bsxfun(@times,et{3}{1},x{1}); %sinF2F1tanX
for m=1:2
    x{m}=ifftdim(x{m},2);       
    x{m}=fftdim(x{m},1);%x{1}=F1F2HsinF2F1tanX
end
x{1}=bsxfun(@times,etg{2}{1},x{1})+bsxfun(@times,et{2}{1},x{2});%tan'sin tan +tan( + )
x{1}=ifftdim(x{1},1);

x{1}=fftdim(x{1},3);
GC{1}=x{1};
x{1}=bsxfun(@times,x{1},et{2}{2});%tan2
x{1}=ifftdim(x{1},3);
x{1}=fftdim(x{1},1);
x{1}=bsxfun(@times,x{1},et{3}{2});
x{1}=ifftdim(x{1},1);
x{1}=fftdim(x{1},3);
x{1}=bsxfun(@times,x{1},et{2}{2});
x{1}=ifftdim(x{1},3);

x{1}=fftdim(x{1},2);
GC{2}=x{1};
x{1}=bsxfun(@times,x{1},et{2}{3});%tan3
x{1}=ifftdim(x{1},2);
x{1}=fftdim(x{1},3);
x{1}=bsxfun(@times,x{1},et{3}{3});
x{1}=ifftdim(x{1},3);
x{1}=fftdim(x{1},2);
x{1}=bsxfun(@times,x{1},et{2}{3});%F2H?
% x{1}=ifftdim(x{1},2);%xc

for m=1:2:3
    x{1}=fftdim(x{1},m);
end
GB{1}=x{1};
x{1}=bsxfun(@times,x{1},et{1});
for m=3:-1:1    
    x{1}=ifftdim(x{1},m);   
end
G{4}=x{1};


%Second rotation
x{1}=bsxfun(@times,xB{2},et{2}{2});%tan
x{2}=bsxfun(@times,xB{2},etg{2}{2});%tan'
for m=1:2
    x{m}=ifftdim(x{m},3); %F3'tan   F3'tan'
    x{m}=fftdim(x{m},1); %F1F3'tan   F1F3'tan'
end
x{2}=bsxfun(@times,etg{3}{2},x{1})+bsxfun(@times,et{3}{2},x{2});%sin'tan +sintan'
x{1}=bsxfun(@times,et{3}{2},x{1});%sintan
for m=1:2
    x{m}=ifftdim(x{m},1);
    x{m}=fftdim(x{m},3);
end
x{1}=bsxfun(@times,etg{2}{2},x{1})+bsxfun(@times,et{2}{2},x{2});
x{1}=ifftdim(x{1},3); %F3'

x{1}=fftdim(x{1},2); %F2
GC{3}=x{1}; %before 3rd rotation
x{1}=bsxfun(@times,x{1},et{2}{3});
x{1}=ifftdim(x{1},2);
x{1}=fftdim(x{1},3);
x{1}=bsxfun(@times,x{1},et{3}{3});
x{1}=ifftdim(x{1},3);
x{1}=fftdim(x{1},2);
x{1}=bsxfun(@times,x{1},et{2}{3});
% x{1}=ifftdim(x{1},2);%xc

for m=1:2:3
    x{1}=fftdim(x{1},m);
end
GB{2}=x{1};
x{1}=bsxfun(@times,x{1},et{1});
for m=3:-1:1    
    x{1}=ifftdim(x{1},m);    
end
G{5}=x{1};


%Third rotation
x{1}=bsxfun(@times,xB{3},et{2}{3}); %tan
x{2}=bsxfun(@times,xB{3},etg{2}{3}); %tan'
for m=1:2
    x{m}=ifftdim(x{m},2); %F2'
    x{m}=fftdim(x{m},3); %F2'F3
end
x{2}=bsxfun(@times,etg{3}{3},x{1})+bsxfun(@times,et{3}{3},x{2}); %sin'tan+sintan'
x{1}=bsxfun(@times,et{3}{3},x{1}); %sintan
for m=1:2
    x{m}=ifftdim(x{m},3); %F3'
    x{m}=fftdim(x{m},2); %F2F3'
end
x{1}=bsxfun(@times,etg{2}{3},x{1})+bsxfun(@times,et{2}{3},x{2});
% x{1}=ifftdim(x{1},2);%xc

for m=1:2:3
    x{1}=fftdim(x{1},m);        
end
GB{3}=x{1};
x{1}=bsxfun(@times,x{1},et{1});
for m=3:-1:1       
    x{1}=ifftdim(x{1},m);        
end
G{6}=x{1};