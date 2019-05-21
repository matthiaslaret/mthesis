function gap_opening

close all;

n=100;

r = exp(linspace(log(0.001*206265*149597870700),log(10*206265*149597870700),n));


h = 0.0001 * r .* (r./149597870700).^(1/3);
alpha = 0.01;
Ms = 3*10^5 * 2e30;
Mb = 40* 2e30;
q = Mb/Ms;
Rh = r * (q/3)^(1/3) ;


fcn = 3/4 .* h./Rh + 50*alpha*h.^2./(q*r.^2) * sqrt(Ms./(Mb + Ms)) ;
fcn = 3/4 .* h./r * (q./3)^(-1/3) + 50*alpha./q * (h./r).^2;

%fcn - fcn1

figure; 
loglog(r/149597870700/206265,fcn,r/149597870700/206265,h./r,r/149597870700/206265,ones(1,n),'--');
hold off
xlabel('radius in disk [pc]')
ylabel('gap-opening paramater g')
grid on;

end

