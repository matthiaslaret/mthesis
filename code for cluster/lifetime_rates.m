function [lifetime]=lifetime_rates(M,mu,m1,q,a0,e0,alpha,sigma,cs)

%constants
G = 6.67408e-11;
c = 299792458;

options = odeset('RelTol',1e-6);

[tv,Yv]=ode23tb(@(t,y) funsys_rates(t,y,G,c,M,mu,m1,q,alpha,sigma,cs),[0 1e20],[a0;e0],options);

lifetime = tv(end)/31536000;

end