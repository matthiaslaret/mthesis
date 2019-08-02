function rates_galaxy(Ir,Im,Iq,Ia,Ie)

%% CONSTANTS
SM = 2e30; %Solar Mass
G = 6.67408e-11;
c = 299792458;
AU = 149597870700; % AU in meters
pc = 206265*AU; %pc in meters


%
%% GALAXY PARAMETERS
Msmbh = 1e8 * SM ;
gamma = 1.5;
beta = 3.2;

Ntotalgal = 1e6;
fb = 0.01; %binary fraction
fd = 0.01; % down to disk fraction

alpha = 0.1;
aspectratio = 0.01;
sigma = ( Msmbh/(3.7e15) )^(1/4.38) ; %velocity dispersion m/s, using sigma-M relation
fg = 0.1; %gas surface density divided by total surface density
rb = G*Msmbh /sigma^2;

SSR = 2*G*Msmbh/c^2;

%
%% MIN MAX PARAMETERS
% RADIUS AGN
Rmin = 100*SSR; %inner radius of AGN
Rmax = 1*pc; %outer radius of AGN

% STELlAR MASS
msmax = 200*SM; %max stellar mass
msmin = 0.01*SM; %min stellar mass
mscr = 5 * SM; % critical stellar mass

% BH MASS
mbhmax = 50 * SM; %max BH mass
mbhmin = 5 * SM; %min BH mass 

% MASS RATIO 
qmin = 0.01;
qmax = 1; 

% Orbital sep
a0min = 0.1 * AU;
a0max = 100 * AU;

% Eccentricity
e0min = 0.001;
e0max = 0.95;

%
%% PARTITIONS

%Ir = 5; %numer of radial bins
%Im = 5; %number of mass bins
%Iq = 5; %number of mass ratio bins
%Ia = 5; %number of orbital separation bins
%Ie = 5; %number of orbital ecentricity bins


%
%% PARAMETERS, NO INPUT

k0 = Ntotalgal/SM/((0.56 - 1/(0.7)*(msmin/SM)^(0.7) - 0.04/(1.3)*(msmax/SM)^(-1.3) ));
Mtot = k0*SM^2*(0.22 - 1/(1.7)*(msmin/SM)^(1.7) - 0.04/(0.3)*(msmax/SM)^(-0.3) );
rho0 = Mtot/(4*pi*rb^3) * (1/(3-gamma) - 1/(3-beta) )^(-1) ;  
k = 1/SM * (- 0.04/(1.3)* ((msmax/SM)^(-1.3) - (mscr/SM)^(-1.3) ))/(0.22 - 1/(1.7)*(msmin/SM)^(1.7) - 0.04/(0.3)*(msmax/SM)^(-0.3) );
xi0 = - 1.35 /(mbhmax^(-1.35) - mbhmin^(-1.35));


delr = (Rmax - Rmin)/Ir;
delm = (mbhmax - mbhmin)/Im ;
delq = (qmax - qmin)/Iq;
dela = (a0max - a0min)/Ia;
dele = (e0max - e0min)/Ie;


%
%% MERGER RATE DENSITY
totalrate=[];
totalmrgtimes=[];

for ll=0:(Ir-1)
    
r1 = Rmin + ll*delr;
r2 = Rmin + (ll+1)*delr; 
r = Rmin + ll*delr + 0.5*delr ; %average radius of [r1,r2]

[ratesshell,mrgtimesshell] = rates_radialshell(r1,r2,r);
totalrate = [totalrate ratesshell];
totalmrgtimes = [totalmrgtimes, mrgtimesshell];

end

MERGERTIME_AVERAGE = mean(totalmrgtimes)
MERGERTIME_MEDIAN = median(totalmrgtimes)
RATE = sum(totalrate)

filename = strcat(datestr(now, 'dd_mm_yy_HH_MM'),'.mat'); 
save(filename)

%
%% MERGER RATE FOR RADIAL SHELL [r1,r2]
function [mergerrate,mrgtimes] = rates_radialshell(r1,r2,r)
       
% BH in [r1,r2]
if  r1 > rb
    Ntot = k/(3-beta)*4*pi*rho0*rb^3*( (r2/rb).^(3-beta) - (r1/rb).^(3-beta) ); 
else
    Ntot = k/(3-gamma)*4*pi*rho0*rb^(gamma)*( (r2)^(3-gamma) - r1^(3-gamma) );
end

NBBH = Ntot*fb*fd/2;

%% PARTITIONS
%primary mass

for i=0:(Im-1)
% number of primary BBH in mass bin i     
if  r1 > rb
    Nm(i+1) = 0.5*fb*fd*4*pi*k*rho0*xi0/1.35*((mbhmin + i*delm)^(-1.35) - (mbhmin +(i+1)*delm)^(-1.35))*rb^3/(3-beta)*( (r2/rb).^(3-beta) - (r1/rb).^(3-beta) );
else
    Nm(i+1) = 0.5*fb*fd*4*pi*k*rho0*xi0/1.35*((mbhmin + i*delm)^(-1.35) - (mbhmin +(i+1)*delm)^(-1.35))*rb^(gamma)*1/(3-gamma)*( (r2)^(3-gamma) - r1^(3-gamma) );
end
% average mass in mass bin i
m(i+1) = 1.35/0.35 * ((mbhmin + i*delm)^(-0.35) - (mbhmin +(i+1)*delm)^(-0.35))/((mbhmin + i*delm)^(-1.35) - (mbhmin +(i+1)*delm)^(-1.35));
end


%mass ratio

for i=0:(Iq-1)
phi0 = NBBH/(qmax - qmin);
Nq(i+1) = phi0 *delq;
q(i+1) = (i+0.5)*delq + qmin;
end

%orbital separation

for i=0:(Ia-1)
chi0 = NBBH/log(a0max/a0min);
Na(i+1) = chi0 * log((a0min + (i+1)*dela )/(a0min + i*dela)) ;
a(i+1) = dela/log((a0min + (i+1)*dela )/(a0min + i*dela));
end

%orbital eccentricity

for i=0:(Ie-1)
zeta0 = 2*NBBH/(e0max^2 - e0min^2);
Ne(i+1) = zeta0 * (i+0.5)*dele^2 + zeta0 * e0min *dele ;
e(i+1) = (6*e0min^2 + 6*e0min*(2*i*dele + dele) + 2*dele^2*(3*i^2 + 3*i +1) )/(3*(2*e0min + 2*i * dele + dele)) ;
end

%
%% LIFETIMES

sigmar = fg * sigma^2./(pi*r*G); %gas surface density at r 
csr = aspectratio * sqrt(G*Msmbh/r) ; %sound speed at r 



mrgtimes = [];
contr = [];
for i=1:Im
    m1 = m(i);
    [r/Rmax,i]
    for j=1:Iq
        q1 = q(j);
        Mbin = m1 + m1*q1;
        mu = Mbin*(q1/(1+q1)^2);
        for l=1:Ia
            a0 = a(l);
%             a0 = 20*149597870700;
            for n=1:Ie
                e0 = e(n);
%                 e0 = 0.001;
                
                mergertime = lifetime_rates(Mbin,mu,m1,q1,a0,e0,alpha,sigmar,csr);
                mrgtimes = [mrgtimes mergertime];
                ratecontri = (1/mergertime)*Nm(i)*Nq(j)*Na(l)*Ne(n);
                contr = [contr ratecontri];
            end
        end
    end
end

mergerrate = sum(contr)/((NBBH)^3);
plot(mrgtimes,'o');

end
end