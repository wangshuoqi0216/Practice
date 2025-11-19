%% AC equation
clear; clc
syms x y epsilon M t

phi = cos(x).*cos(y).*exp(-t);

lap_phi = diff(phi,x,2) + diff(phi,y,2);

fun = (phi.^3-phi);

mu = - epsilon^2 * lap_phi + fun;

rhs = simplify(diff(phi,t) + M*mu);

rhs=char(rhs);
rhs= strrep(rhs,'*','.*');
rhs= strrep(rhs,'/','./');
rhs= strrep(rhs,'^','.^')

