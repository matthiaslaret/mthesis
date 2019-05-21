function orbit_evolution

%parameters binary
M = 14*2e30; % total binary mass in solarmasses
mu = M/4; % reduced mass (mass ratio = 1)

a0 = 20*149597870700; % in AU
e0 = 0.3;

%Hulse-Taylor
%--------------
% a0 = 0.0130 * 149597870700;
% e0 = 0.6171334;
% M = 2.8*2e30;

%parameters disk
R = 100*a0;   % disk radius
H = 0.1 * R; % disk height
H = 10*149597870700; % in AU
alpha = 1*1e-2;   % viscosity parameter 
Md = 50*2e30; % disk  mass in solarmasses
sigma = 2e+04; % surface density kg/m^2
cs = 0.05 * 0.5 * 200e3; %sound speed

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

% % as a function of eccentricity
% =======================================
% =======================================
e0=linspace(0.01,0.98,2);

dec=0;
x=0;

for i=1:length(e0)
lifetimes_0(i) = lifetime(x,dec,M,mu,a0,e0(i),R,H,alpha,Md,sigma,cs);
end

x=1;
for i=1:length(e0)
lifetimes_gas(i) = lifetime(x,dec,M,mu,a0,e0(i),R,H,alpha,Md,sigma,cs);
end

dec=1;
for i=1:length(e0)
timeofdec(i) = lifetime(x,dec,M,mu,a0,e0(i),R,H,alpha,Md,sigma,cs);
end

dec=2;
for i=1:length(e0)
eccatdec(i) = lifetime(x,dec,M,mu,a0,e0(i),R,H,alpha,Md,sigma,cs);
end

dec=3;
for i=1:length(e0)
aatdec(i) = lifetime(x,dec,M,mu,a0,e0(i),R,H,alpha,Md,sigma,cs);
end


%comparison GW vs GW+disk
semilogy(e0,lifetimes_0,e0,lifetimes_gas,e0,lifetimes_gas - timeofdec,'Marker','.','MarkerSize',20,'Linewidth',3);
grid
legend('GW-driven merger time', '(GW + Disk)-driven merger time','Merger time from Decoupling')
xlabel('Initial Eccentricity','FontSize',18,'FontWeight','bold')
ylabel('Time [years]','FontSize',18,'FontWeight','bold');
txt = ['Binary Mass = ',num2str(M/2e30),' SM; \mu = ',num2str(mu/2e30),' SM; Initial a = ',num2str(a0/149597870700),' AU; sound speed c_s = ',num2str(cs/1000),' km/s; CBD \Sigma = ',num2str(sigma),' kg/m^2; \alpha = ',num2str(alpha)];
title(txt)


%decoupling a & e
[ax,h1,h2]=plotyy(e0,aatdec/149597870700,e0,eccatdec);
grid
xlabel('Initial Eccentricity','FontSize',18,'FontWeight','bold')
ylabel(ax(1), 'Semi-Major Axis at Decoupling [AU]','FontSize',18,'FontWeight','bold');
ylabel(ax(2), 'Eccentricity at Decoupling','FontSize',18,'FontWeight','bold');
set(h1,'Marker','.','MarkerSize',20,'Linewidth',3)
set(h2,'Marker','.','MarkerSize',20,'Linewidth',3)
set(ax,'fontsize',15,'FontWeight','bold','LineWidth',1.5)
txt = ['Binary Mass = ',num2str(M/2e30),' SM; \mu = ',num2str(mu/2e30),' SM; Initial a = ',num2str(a0/149597870700),' AU; CBD scale height = ',num2str(H/149597870700),' AU; CBD \Sigma = ',num2str(sigma),' kg/m^2; \alpha = ',num2str(alpha)];
title(txt)

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

% % as a function of fixed mass ratio Md/Mb with variable M/mu/Md
% % =======================================
% % =======================================
% M = linspace(2*2e30,200*2e30,3);
% mu = M/4;
% Md = M/200;
% 
% 
% dec=0;
% x=0;
% 
% for i=1:length(M)
% lifetimes_0(i) = lifetime(x,dec,M(i),mu(i),a0,e0,R,H,alpha,Md(i));
% end
% 
% x=1;
% for i=1:length(M)
% lifetimes_gas(i) = lifetime(x,dec,M(i),mu(i),a0,e0,R,H,alpha,Md(i));
% end
% 
% dec=1;
% for i=1:length(M)
% timeofdec(i) = lifetime(x,dec,M(i),mu(i),a0,e0,R,H,alpha,Md(i));
% end
% 
% dec=2;
% for i=1:length(M)
% eccatdec(i) = lifetime(x,dec,M(i),mu(i),a0,e0,R,H,alpha,Md(i));
% end
% 
% dec=3;
% for i=1:length(M)
% aatdec(i) = lifetime(x,dec,M(i),mu(i),a0,e0,R,H,alpha,Md(i));
% end

% %comparison GW vs GW+disk
% [ax,h1,h2]=plotyy(M/2e30,lifetimes_0,M/2e30,lifetimes_gas,'semilogy','semilogy');
% grid
% legend([h1 h2],'GW-driven', '(GW + Disk)-driven')
% xlabel('Binary Mass [SM]','FontSize',18,'FontWeight','bold')
% ylabel(ax(1), 'Merger Time [years]','FontSize',18,'FontWeight','bold');
% ylabel(ax(2), 'Merger Time [years]','FontSize',18,'FontWeight','bold');
% set(h1,'Marker','.','MarkerSize',20,'Linewidth',3)
% set(h2,'Marker','.','MarkerSize',20,'Linewidth',3)
% set(ax,'fontsize',15,'FontWeight','bold','LineWidth',1.5)
% txt = ['Initial a = ',num2str(a0/149597870700),' AU; CBD Radius = ',num2str(R/149597870700),' AU; CBD Height = ',num2str(H/149597870700),' AU; \alpha = ',num2str(alpha)];
% title(txt)
% 
% %decoupling
% [ax,h1,h2]=plotyy(e0,lifetimes_gas - timeofdec,e0,eccatdec);
% grid
% xlabel('Binary Mass','FontSize',18,'FontWeight','bold')
% ylabel(ax(1), '(Merger Time - Decoupling Time) [years]','FontSize',18,'FontWeight','bold');
% ylabel(ax(2), 'Eccentricity at Decoupling','FontSize',18,'FontWeight','bold');
% set(h1,'Marker','.','MarkerSize',20,'Linewidth',3)
% set(h2,'Marker','.','MarkerSize',20,'Linewidth',3)
% set(ax,'fontsize',15,'FontWeight','bold','LineWidth',1.5)
% txt = ['Binary Mass = ',num2str(M/2e30),' SM; \mu = ',num2str(mu/2e30),' SM; Initial a = ',num2str(a0/149597870700),' AU; CBD Radius = ',num2str(R/149597870700),' AU; CBD Height = ',num2str(H/149597870700),' AU; CBD Mass = ',num2str(Md/2e30),' SM; \alpha = ',num2str(alpha)];
% title(txt)
% 
% %decoupling a & e
% [ax,h1,h2]=plotyy(e0,aatdec,e0,eccatdec);
% grid
% xlabel('Initial Eccentricity','FontSize',18,'FontWeight','bold')
% ylabel(ax(1), 'Semi-Major Axis at Decoupling','FontSize',18,'FontWeight','bold');
% ylabel(ax(2), 'Eccentricity at Decoupling','FontSize',18,'FontWeight','bold');
% set(h1,'Marker','.','MarkerSize',20,'Linewidth',3)
% set(h2,'Marker','.','MarkerSize',20,'Linewidth',3)
% set(ax,'fontsize',15,'FontWeight','bold','LineWidth',1.5)
% txt = ['Binary Mass = ',num2str(M/2e30),' SM; \mu = ',num2str(mu/2e30),' SM; Initial a = ',num2str(a0/149597870700),' AU; CBD Radius = ',num2str(R/149597870700),' AU; CBD Height = ',num2str(H/149597870700),' AU; CBD Mass = ',num2str(Md/2e30),' SM; \alpha = ',num2str(alpha)];
% title(txt)


end

