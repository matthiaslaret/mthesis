function stoneplot_v2

% find a_0 at given r s.t. mergerime is less than t_want

close all;

t_want_array = [1e6 1e8];
for k=1:length(t_want_array)

t_want = t_want_array(k); %lifetime in years

%general parameters
G = 6.6740e-11;

e0 = 0.8;
M = 10*2e30; % total binary mass in solarmasses
mu = M/4; % reduced mass (mass ratio = 1)


n=2;
r = exp(linspace(log(0.01*206265*149597870700),log(10*206265*149597870700),n)); 


for i=1:n

a0_lin=exp(linspace(log(0.05*149597870700),log(50*149597870700),10));
for j=1:length(a0_lin)
    a0=a0_lin(j);

%disk parameters TO
Ms = 3*10^6 * 2e30; %central mass
h = 2^(-3/2) * r * 0.05; %scale height
sigma = (200e3)^2./(pi*r*G); %surface density
alpha = 0.01;
Md =0; %only needed for lifetimes.m to work
cs = 0.05 * 0.5 * 200e3; %sound speed

dec=0; % no decoupling
% x=0; %no gas
% lifetimes_0(j) = lifetime(x,dec,M,mu,a0,e0,r(i),h(i),alpha,Md,sigma(i),cs);

x=1; %with gas
lifetimes_1(j) = lifetime(x,dec,M,mu,a0,e0,r(i),h(i),alpha,Md,sigma(i),cs);

end

%search gas
ind = find((lifetimes_1 - t_want) < 0);
a0_want(i) = a0_lin(min(ind))/149597870700;


end

a0_array(k,:) = a0_want;

end

loglog(r/149597870700/206265,a0_array(1,:),'--o',r/149597870700/206265,a0_array(2,:),'--o','LineWidth',3)
xlabel('radius [pc]');
ylabel('Semi-major axis a_0 [AU]')
txt = ['Binary Mass = ',num2str(M/2e30),' SM; \mu = ',num2str(mu/2e30),' SM; Initial e = ',num2str(e0),'; sound speed c_s = ',num2str(cs/1000),' km/s; \alpha = ',num2str(alpha)];
title(txt)
grid on;



end