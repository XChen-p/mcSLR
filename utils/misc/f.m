function varargout=f(varargin)
    if nargin == 0
        if nargout == 1
            varargout{1}  =   figure();
        else
            figure();
            datacursormode on;
        end
    else
        if ~isnumeric(varargin{1})
            varargin{1}=str2num(varargin{1});
        end
        figure(varargin{1});
        datacursormode on;
    end
    set(gcf, 'InvertHardCopy', 'off');
    drawnow;commandwindow;
end
