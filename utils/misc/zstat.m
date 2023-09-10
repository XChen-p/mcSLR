function [zt, zf] = zstat(varargin)

[tt, dof, ff, dof1, dof2]   =   tstat(varargin{:});

zt  =   t_to_z(tt, dof);
zf  =   f_to_z(ff, dof1, dof2);
