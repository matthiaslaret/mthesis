function  rates_distributions

SM = 2e30;
G = 6.67408e-11;
Msmbh = 1e5 * SM ;
gamma = 1.5;
beta = 3.2;
sigma = 30e3; %m/s
msmax = 50000*SM;
msmin = 0.01*SM;
mscr = 0.05 * SM;
mbhmax = 50 * SM;
mbhmin = 1*SM;
rb = G*Msmbh /sigma^2;


k = 1/SM * (0.56 - 1/(0.7)*(mscr/SM)^(0.7) - 0.04/(1.3)*(msmax/SM)^(-1.3) )/(0.22 - 1/(1.7)*(msmin/SM)^(1.7) - 0.04/(0.3)*(msmax/SM)^(-0.3) );
xi0 = - 1.35 /(mbhmax^(-1.35) - mbhmin^(-1.35));



n=1000;
r = exp(linspace(log(0.01*206265*149597870700),log(10*206265*149597870700),n));
dr = 0.0001*206265*149597870700;


m = 5*SM;

for i=1:n

if r(i) > rb
    dNdm(i) = k.*fBH(m).*2*Msmbh*((3-gamma)/(3-beta))*( ((r(i)+dr)/rb).^(3-beta) - (r(i)/rb).^(3-beta) );
    dN(i) = k*2*Msmbh*((3-gamma)/(3-beta))*( ((r(i)+dr)/rb).^(3-beta) - (r(i)/rb).^(3-beta) ); 
else
    dNdm(i) = k.*fBH(m).*2*Msmbh*rb^(gamma -3)*( (r(i)+dr)^(3-gamma) - r(i)^(3-gamma) );
    dN(i) = k*2*Msmbh*rb^(gamma -3)*( (r(i)+dr)^(3-gamma) - r(i)^(3-gamma) );
end

end

m = 50*SM;

for i=1:n

if r(i) > rb
    dNdm2(i) = k.*fBH(m).*2*Msmbh*((3-gamma)/(3-beta))*( ((r(i)+dr)/rb).^(3-beta) - (r(i)/rb).^(3-beta) );
else
    dNdm2(i) = k.*fBH(m).*2*Msmbh*rb^(gamma -3)*( (r(i)+dr)^(3-gamma) - r(i)^(3-gamma) );
end

end

loglog(r/206265/149597870700,dNdm,r/206265/149597870700,dN,r/206265/149597870700,dNdm2,'--')
xlabel('radius [pc]');
ylabel('dN/dm');
grid on;


    function [fBH] = fBH(m)
        fBH = xi0*m.^(-2.35);
    end


end

