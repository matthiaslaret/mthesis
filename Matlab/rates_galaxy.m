function rates_galaxy


SM = 2e30;
G = 6.67408e-11;
Msmbh = 1e6 * SM ;
gamma = 1.5;
beta = 3.2;
sigma = 30e3; %m/s
msmax = 50000*SM;
msmin = 0.01*SM;
mscr = 5 * SM;
mbhmax = 50 * SM;
mbhmin = 5*SM;
rb = G*Msmbh /sigma^2;


Ntotalgal = 2e2;
k0 = Ntotalgal/SM/((0.56 - 1/(0.7)*(mscr/SM)^(0.7) - 0.04/(1.3)*(msmax/SM)^(-1.3) ));
Mtot = k0*SM^2*(0.22 - 1/(1.7)*(msmin/SM)^(1.7) - 0.04/(0.3)*(msmax/SM)^(-0.3) );
rho0 = Mtot/(4*pi*rb^3) * (1/(3-gamma) - 1/(3-beta) )^(-1) ;  
k = 1/SM * (0.56 - 1/(0.7)*(mscr/SM)^(0.7) - 0.04/(1.3)*(msmax/SM)^(-1.3) )/(0.22 - 1/(1.7)*(msmin/SM)^(1.7) - 0.04/(0.3)*(msmax/SM)^(-0.3) );
xi0 = - 1.35 /(mbhmax^(-1.35) - mbhmin^(-1.35));

%radial partition
Ir = 2; %numer of radial bins
delr = 10*206265*149597870700/Ir;

totalrate=[];

for ll=0:(Ir-1)
    
r1 = ll*delr;
r2 = (ll+1)*delr; 
dr = r2-r1;
ratesshell = rates_radialshell(r1,dr)
totalrate = [totalrate ratesshell];

end

sum(totalrate)

function mergerrate = rates_radialshell(r,dr)
% merger rate in a radial shell [r,r+dr]

% r = 0.01*206265*149597870700;
% dr = 0.0001*206265*149597870700;

%primary mass
Im = 2; %number of mass bins
delm = (mbhmax - mbhmin)/Im ;

for i=0:(Im-1)
if  r > rb
    Nm(i+1) = 0.5*4*pi*k*rho0*xi0/1.35*((mbhmin + i*delm)^(-1.35) - (mbhmin +(i+1)*delm)^(-1.35))*rb^3/(3-beta)*( ((r+dr)/rb).^(3-beta) - (r/rb).^(3-beta) );
    Ntot = k/(3-beta)*4*pi*rho0*rb^3*( ((r+dr)/rb).^(3-beta) - (r/rb).^(3-beta) ); 
else
    Nm(i+1) = 0.5*4*pi*k*rho0*xi0/1.35*((mbhmin + i*delm)^(-1.35) - (mbhmin +(i+1)*delm)^(-1.35))*rb^(gamma)*1/(3-gamma)*( (r+dr)^(3-gamma) - r^(3-gamma) );
    Ntot = k/(3-gamma)*4*pi*rho0*rb^(gamma)*( (r+dr)^(3-gamma) - r^(3-gamma) );
end
m(i+1) = 1.35/0.35 * ((mbhmin + i*delm)^(-0.35) - (mbhmin +(i+1)*delm)^(-0.35))/((mbhmin + i*delm)^(-1.35) - (mbhmin +(i+1)*delm)^(-1.35));
end


%mass ratio
Iq = 2; %number of mass ratio bins
delq = 0.9/Iq;

for i=0:(Iq-1)
Nq(i+1) = 10/18 *delq *Ntot ;
q(i+1) = (i+0.5)*delq + 0.1;
end

%orbital separation
Ia = 2; %number of orbital separation bins
dela = 499.99*149597870700/Ia;

for i=0:(Ia-1)
Na(i+1) = (Ntot/(2*log(50000)))* log((0.01*149597870700 + (i+1)*dela )/(0.01*149597870700 + i*dela)) ;
a(i+1) = dela/log((0.01*149597870700 + (i+1)*dela )/(0.01*149597870700 + i*dela));
end

%orbital eccentricity
Ie = 2; %number of orbital ecentricity bins
dele = 1/Ie;

for i=0:(Ie-1)
Ne(i+1) = Ntot*(i+0.5)*dele^2 ;
e(i+1) = (((i+1)*dele)^3 - (i*dele)^3)/(3*(i+0.5)*dele^2) ;
end


%lifetimes
x=1;
dec=0;
Mbin = m(1) + q(1)*m(1);
mu = Mbin*(q(1)/(1+q(1))^2);
a0 = a(1);
e0 = e(1);
R = r; %not needed
H = r; %not needed
alpha = 0.01;
Md = r; %not needed
r = r1 + 0.5*dr; %average radius of mass bin for disk property
sigma = (200e3)^2./(pi*r*G); %surface density
cs = 0.05 * 0.5 * 200e3; %sound speed

mrgtimes = [];
contr = [];
for i=1:Im
    for j=1:Iq
        Mbin = m(i) + m(i)*q(j);
        mu = Mbin*(q(j)/(1+q(j))^2);
        for l=1:Ia
            a0 = a(l);
            for n=1:Ie
                [i,j,l,n]
                e0 = e(n);
                mergertime = lifetime(x,dec,Mbin,mu,a0,e0,R,H,alpha,Md,sigma,cs);
                mrgtimes = [mrgtimes mergertime];
                ratecontri = (1/mergertime)*Nm(i)*Nq(j)*Na(l)*Ne(n);
                contr = [contr ratecontri];
            end
        end
    end
end

mergerrate = sum(contr)/((0.5*Ntot)^3)
plot(mrgtimes,'o')

end

end