%   Replacement for plot that automatically plots real and imag, and squeezes
function plt(varargin)

h = ishold;

if nargin == 1
    builtin('plot',squeeze(real(varargin{1})));
    if ~isreal(varargin{1})
        hold on;
        builtin('plot',squeeze(imag(varargin{1})));
    end
elseif nargin == 2
    if isnumeric(varargin{2})
        builtin('plot',varargin{1}, squeeze(real(varargin{2})));
        if ~isreal(varargin{2})
            hold on;
            builtin('plot',varargin{1}, squeeze(imag(varargin{2})));
        end
    else
        plot(1:length(varargin{1}), squeeze(real(varargin{1})), varargin{2});
    end
else
    if isnumeric(varargin{2})
        builtin('plot',varargin{1}, squeeze(real(varargin{2})), varargin{3:end});
        if ~isreal(varargin{2})
            hold on;
            builtin('plot',varargin{1}, squeeze(imag(varargin{2})), varargin{3:end});
        end
    else
        plot(1:length(varargin{1}), squeeze(real(varargin{1})), varargin{2:end});
    end
end
if ~h
    hold off;
end
drawnow;
commandwindow;
