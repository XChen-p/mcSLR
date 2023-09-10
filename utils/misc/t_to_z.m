function z = t_to_z(t, dof)

    z   =   -sign(t).*norminv(tcdf(abs(t), dof, 'upper'));
