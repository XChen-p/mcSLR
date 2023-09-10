%   Replacement for imshow that automatically squeezes input image
function show(img, varargin)

if size(img,1)==numel(img)
    if mod(size(img,1), sqrt(size(img,1))) == 0
        img =   reshape(img, sqrt(size(img,1)), []);
    end
end

datacursormode on;
if nargin == 1
    imshow(squeeze(img), [], 'colormap', jet);
else
    imshow(squeeze(img), varargin{:});
end
%set(gca,'YDir','Normal');
drawnow;
commandwindow;
