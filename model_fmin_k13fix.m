function [ss,sol,rout] = model_fmin_k13fix(pars)

global data

xdc = data.x;
rdc = data.r;

k2   = pars(1);
rin  = rdc(1);
rout = rdc(end);

cst  =  exp(-k2*xdc(end));
k3   = (rout-rin*cst)/(1-cst);
k1   = rin-k3;
% k1 + k3 = rin <=> k1 = rin - k3
% rout = k1*c2 + k3 <=> rout = (rin-k3)*cst + k3 = rin*cst + k3 (1-cst)
% k3 = (rout-rin*cst)/(1-cst)
% k1 = rin - k3;
rad = k1*exp(-k2*xdc)+k3;

rout = (rad - rdc);
ss = rout'*rout;
