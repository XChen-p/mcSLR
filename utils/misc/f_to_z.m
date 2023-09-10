function z = f_to_z(ff, dof1, dof2);

    z   =   -sign(ff).*norminv(fcdf(abs(ff), dof1, dof2, 'upper'));
