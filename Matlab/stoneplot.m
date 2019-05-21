function stoneplot

%%% means a0 vs mergertime

close all;

%general parameters
G = 6.6740e-11;
n=2;
r = exp(linspace(log(0.01*206265*149597870700),log(10*206265*149597870700),n)); 


%binary parameters
a0 = 1*149597870700; % in AU

% % % a0_lin=exp(linspace(log(0.1*149597870700),log(10*149597870700),20));
% % % for j=1:length(a0_lin)
% % %     a0=a0_lin(j);

e0 = 0.8;
M = 10*2e30; % total binary mass in solarmasses
mu = M/4; % reduced mass (mass ratio = 1)


%disk parameters TO
Ms = 3*10^6 * 2e30; %central mass
h = 2^(-3/2) * r * 0.05; %scale height
sigma = (200e3)^2./(pi*r*G); %surface density
alpha = 0.01;
Md =0; %only needed for lifetimes.m to work
cs = 0.05 * 0.5 * 200e3; %sound speed

dec=0; % no decoupling
x=0; %no gas
for i=1:n
    lifetimes_0(i) = lifetime(x,dec,M,mu,a0,e0,r(i),h(i),alpha,Md,sigma(i),cs);
end

x=1; %with gas
for i=1:n
    lifetimes_1(i) = lifetime(x,dec,M,mu,a0,e0,r(i),h(i),alpha,Md,sigma(i),cs);
end

loglog(r/149597870700/206265,lifetimes_0,'--o',r/149597870700/206265,lifetimes_1,'--o','LineWidth',3)
xlabel('radius [pc]')
ylabel('merger time [years]')
legend('GW-driven','(GW + disk)-driven')
txt = ['Binary Mass = ',num2str(M/2e30),' SM; \mu = ',num2str(mu/2e30),' SM; Initial a = ',num2str(a0/149597870700),' AU; Initial e = ',num2str(e0),'; sound speed c_s = ',num2str(cs/1000),' km/s; \alpha = ',num2str(alpha)];
title(txt)
grid on;

% % % minima(j) = min(lifetimes_1);
% % % maxima(j) = max(lifetimes_1);
% % % end
% % % loglog(a0_lin/149597870700,minima,'--o',a0_lin/149597870700,maxima,'--o','LineWidth',3)
% % % xlabel('Initial semi-major axis a_0 [AU]');
% % % ylabel('Merger Time [years]')
% % % legend('at 0.01 pc','at 10 pc')
% % % txt = ['Binary Mass = ',num2str(M/2e30),' SM; \mu = ',num2str(mu/2e30),' SM; Initial e = ',num2str(e0),'; sound speed c_s = ',num2str(cs/1000),' km/s; \alpha = ',num2str(alpha)];
% % % title(txt)
% % % grid on;

figure;
loglog(r/149597870700/206265,sigma,r/149597870700/206265,cs/1e3*ones(1,n),'LineWidth',3)
xlabel('radius [pc]')
legend('Surface density \Sigma [kg/m^2]','soundspeed [km/s]')
grid on;

end