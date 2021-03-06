function Fv = funsys_3(t,Y,G,c,M,mu,q,m1,sigma)
% Y(1) = a 
% Y(2) = e

% GW driven
Fv(1,1) = - 10^12 * G^3/c^5 * 64/5 * (mu*M^2)/Y(1)^3 * 1/(1-Y(2)^2)^(7/2) *( 1 + 73/24 * Y(2)^2 + 37/96 * Y(2)^4);
Fv(2,1) = - 10^12 * G^3/c^5 * 304/15 * (mu*M^2)/Y(1)^4 * Y(2)/(1-Y(2)^2)^(5/2) *( 1 + 121/304 * Y(2)^2);

% GW plus disk

bracket1 = (Y(1)-3*(Y(1)-Y(1)/3))/(6*((Y(1)-Y(1)/3) - Y(1) )^3) ;
bracket2 = (Y(1))/(-6*(Y(1))^3) ;
bracket3 = (Y(1)-3*(10000*Y(1)))/(6*(10000*Y(1) - Y(1) )^3) ;
bracket4 = (Y(1)-3*(Y(1)+Y(1)/3))/(6*((Y(1)+Y(1)/3) - Y(1) )^3) ;

bracket5 = (Y(1)-4*(Y(1)-Y(1)/3))/(12*(Y(1)-(Y(1)-Y(1)/3))^4) ;
bracket6 = (Y(1))/(12*(Y(1))^4) ;
bracket7 = (Y(1)-4*(10000*Y(1)))/(12*(Y(1) - 10000*Y(1))^4) ;
bracket8 = (Y(1)-4*(Y(1)+Y(1)/3))/(12*(Y(1)-(Y(1)+Y(1)/3)  )^4) ;


Fv1 = 2*pi*0.8*q/m1 * sqrt(G*m1) * Y(1)^(7/2) *sigma*(bracket1 - bracket2 - bracket3 + bracket4 );
Fv2 = 2*pi*0.07*q/m1 * sqrt(G*m1) * Y(1)^(7/2)*Y(2) *sigma*(bracket5 - bracket6 - bracket7 + bracket8 );

% 
% Fv1 = 2*5.6*q/m1*sigma*Y(1)^2*sqrt(G*m1)*Y(1)*Y(1)^(-3/2)*2^3 ; 
% Fv2 = 2*0.59*q/m1*sigma*Y(1)^2*sqrt(G*m1)*Y(1)^(-3/2)*2^4*Y(2) ;
    
Fv(1,1) =  - G^3/c^5 * 64/5 * (mu*M^2)/Y(1)^3 * 1/(1-Y(2)^2)^(7/2) *( 1 + 73/24 * Y(2)^2 + 37/96 * Y(2)^4) + Fv1;
Fv(2,1) =  - G^3/c^5 * 304/15 * (mu*M^2)/Y(1)^4 * Y(2)/(1-Y(2)^2)^(5/2) *( 1 + 121/304 * Y(2)^2) + Fv2;


% if -Fv1 >  G^3/c^5 * 64/5 * (mu*M^2)/Y(1)^3 * 1/(1-Y(2)^2)^(7/2) *( 1 + 73/24 * Y(2)^2 + 37/96 * Y(2)^4)
%     Fv(1,1) = Fv1; 
%     Fv(2,1) = Fv2;
% end
% 
% if -Fv1 <  G^3/c^5 * 64/5 * (mu*M^2)/Y(1)^3 * 1/(1-Y(2)^2)^(7/2) *( 1 + 73/24 * Y(2)^2 + 37/96 * Y(2)^4)
%     Fv(1,1) = - G^3/c^5 * 64/5 * (mu*M^2)/Y(1)^3 * 1/(1-Y(2)^2)^(7/2) *( 1 + 73/24 * Y(2)^2 + 37/96 * Y(2)^4); 
%     Fv(2,1) = - G^3/c^5 * 304/15 * (mu*M^2)/Y(1)^4 * Y(2)/(1-Y(2)^2)^(5/2) *( 1 + 121/304 * Y(2)^2);
% end

end


