function [out, varargout] = cg(data,st,wi,x0,maxiter, tolgrad, verb)
%
%   out = cg(data, st, wi, x0, maxiter, tolgrad)
%
%   Mark Chiew
%   March 2015
%
%   Uses cg-iterations to solve the inverse NUFFT problem
%   Input:
%           data    -   time-series of non-uniform k-space samples [N, t]
%           st      -   NUFFT structure(s) [t]
%           wi      -   density compensation weights [N, t]
%           x0      -   initial guess
%           maxiter -   maximum iterations
%           tolgrad -   minimum gradient tolerance

nt      =   length(st);
a       =   0.05;
b       =   0.1;

if nargin < 6
    tolgrad =   1E-4;
end
if nargin < 7
    verb = 1;
end

for i = 1:nt
    if verb
        disp(sprintf('Time Point %i',i));
    end
    errnorm =   norm(data(:,i));

    k   =   0;
    x   =   x0(:,:,i);
    g   =   grad(data(:,i), st(i), wi(:,i), x);
    dx  =   -g;

    while norm(g(:)) > tolgrad*errnorm && k < maxiter
        %linesearch
        t   =   1;
        while f(data(:,i),st(i), x+t*dx) > f(data(:,i),st(i),x)+a*t*real(reshape(g,[],1)'*reshape(dx,[],1))
            t   =   b*t;
        end
        x   =   x+t*dx;

        g2  =   grad(data(:,i), st(i), wi(:,i), x);
        dx  =   -g2+dx*(norm(g2(:))/norm(g(:)))^2;
        k   =   k+1;
        g   =   g2;
        if verb
            disp(sprintf('Iter: %03i, Grad: %f, Err: %f',k,norm(g(:))/errnorm,sqrt(f(data(:,i),st(i),x))/errnorm));
        end
    end
    out(:,:,i)  =   x;
    err(i)      =   sqrt(f(data(:,i),st(i),x))/errnorm;
end

if nargout > 0
    varargout{1}    =   err;
end

end

function o = f(data,st,x)
    o   =   norm(reshape(data-nufft(x,st),[],1))^2;
end

function g = grad(data,st,wi,x)
        g   =   2*nufft_adj((nufft(x,st)-data).*wi,st)/(2*numel(data));
end
