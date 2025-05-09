function [ss,sol,rout] = model_fmin(pars)

global data

xdc = data.x;
rdc = data.r;

k1 = pars(1);
k2 = pars(2);
k3 = pars(3);

rad = k1*exp(-k2*xdc)+k3;

rout = (rad - rdc);
ss = rout'*rout;


