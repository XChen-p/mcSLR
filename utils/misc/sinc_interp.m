function y = sinc_interp(t, x, ti)

[TI,T]  =   ndgrid(ti,t);
y       =   sinc(TI-T)*(x-mean(x,1))+mean(x,1);
