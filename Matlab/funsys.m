function Fv = funsys(t,Y,x,G,c,M,mu,H,R,alpha,Md,sigma,cs)
% Y(1) = a 
% Y(2) = e

% GW driven
Fv(1,1) = - 10^12 * G^3/c^5 * 64/5 * (mu*M^2)/Y(1)^3 * 1/(1-Y(2)^2)^(7/2) *( 1 + 73/24 * Y(2)^2 + 37/96 * Y(2)^4);
Fv(2,1) = - 10^12 * G^3/c^5 * 304/15 * (mu*M^2)/Y(1)^4 * Y(2)/(1-Y(2)^2)^(5/2) *( 1 + 121/304 * Y(2)^2);


% GW plus disk
if x==1

Fv(1,1) =  - G^3/c^5 * 64/5 * (mu*M^2)/Y(1)^3 * 1/(1-Y(2)^2)^(7/2) *( 1 + 73/24 * Y(2)^2 + 37/96 * Y(2)^4) -  24*pi*alpha*Y(1)*sigma*cs^2 / (mu*sqrt(G*M/(Y(1)^3)));
Fv(2,1) =  - G^3/c^5 * 304/15 * (mu*M^2)/Y(1)^4 * Y(2)/(1-Y(2)^2)^(5/2) *( 1 + 121/304 * Y(2)^2) + (1/sqrt(1-Y(2)^2) - 1)*12*pi*alpha*cs^2*(1-Y(2)^2)*sigma / (mu*Y(2)*sqrt(G*M/(Y(1)^3))) ;

end

% only disk
if x==2

Fv(1,1) =  - 24*pi*alpha*Y(1)*sigma*cs^2 / (mu*sqrt(G*M/(Y(1)^3)));
Fv(2,1) =  (1/sqrt(1-Y(2)^2) - 1)*12*pi*alpha*cs^2*(1-Y(2)^2)*sigma / (mu*Y(2)*sqrt(G*M/(Y(1)^3)));

end


end


