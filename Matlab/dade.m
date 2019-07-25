function dade 

close all;

SM = 2e30;
G = 6.67408e-11;
c = 299792458; 
alpha = 0.01;
r=1*206265*149597870700;
sigma = (200e3)^2./(pi*r*G); %surface density
cs = 0.05 * 0.5 * 200e3; %sound speed

M = 50*SM;
mu = M/4;
m1=0;
q=0;

%% e
a0 = 10*149597870700;
e = linspace(0.01,0.9,5)

for i=1:5
e0 = e(i);

options = odeset('RelTol',1e-6,'Stats','on');
x=1;


[tv,Yv{i}]=ode23tb(@(t,y) funsys_2(t,y,x,G,c,M,mu,m1,q,alpha,sigma,cs),[0 1e20],[a0;e0],options);
Yv{i}(:,1)=Yv{i}(:,1)/149597870700;


end

semilogy(Yv{1}(:,2),Yv{1}(:,1),Yv{2}(:,2),Yv{2}(:,1),Yv{3}(:,2),Yv{3}(:,1),Yv{4}(:,2),Yv{4}(:,1),Yv{5}(:,2),Yv{5}(:,1))
xlabel('Eccentricity')
ylabel('Semi-Major Axis [AU]')
legendCell = cellstr(num2str(e'));
hleg = legend(legendCell,'Location','southeast','Box','on');
htitle = get(hleg,'Title');
set(htitle,'String','Initial Eccentricity e_0')
grid on;
axis([0 1 10^(-5) 10])

%% a
e0 = 0.3;
a = linspace(1*149597870700,20*149597870700,5);

for i=1:5
a0 = a(i);

options = odeset('RelTol',1e-6,'Stats','on');
x=1;


[tv,Yv{i}]=ode23tb(@(t,y) funsys_2(t,y,x,G,c,M,mu,m1,q,alpha,sigma,cs),[0 1e20],[a0;e0],options);
Yv{i}(:,1)=Yv{i}(:,1)/149597870700;


end

semilogy(Yv{1}(:,2),Yv{1}(:,1),Yv{2}(:,2),Yv{2}(:,1),Yv{3}(:,2),Yv{3}(:,1),Yv{4}(:,2),Yv{4}(:,1),Yv{5}(:,2),Yv{5}(:,1))
xlabel('Eccentricity')
ylabel('Semi-Major Axis [AU]')
legendCell = cellstr(num2str((a/149597870700)'));
hleg = legend(legendCell,'Location','southeast');
htitle = get(hleg,'Title');
set(htitle,'String','Initial Semi-Major Axis a_0 [AU]')
grid on;
axis([0 1 10^(-5) 20])

%% M
e0 = 0.3;
a0 = 10*149597870700;
M = linspace(10,100,5);

for i=1:5
M1 = M(i)*SM;
mu = M1/4;

options = odeset('RelTol',1e-6,'Stats','on');
x=1;


[tv,Yv{i}]=ode23tb(@(t,y) funsys_2(t,y,x,G,c,M1,mu,m1,q,alpha,sigma,cs),[0 1e20],[a0;e0],options);
Yv{i}(:,1)=Yv{i}(:,1)/149597870700;


end

semilogy(Yv{1}(:,2),Yv{1}(:,1),Yv{2}(:,2),Yv{2}(:,1),Yv{3}(:,2),Yv{3}(:,1),Yv{4}(:,2),Yv{4}(:,1),Yv{5}(:,2),Yv{5}(:,1))
xlabel('Eccentricity')
ylabel('Semi-Major Axis [AU]')
legendCell = cellstr(num2str(M'));
hleg = legend(legendCell,'Location','southeast');
htitle = get(hleg,'Title');
set(htitle,'String','Binary Mass M [SM]')
grid on;
axis([0 1 10^(-5) 10])

%% q
e0 = 0.3;
a0 = 10*149597870700;
M = 50*SM;
q = linspace(0.05,1,5);

for i=1:5
q1=q(i);
mu = M*(q1/(1+q1)^2);

options = odeset('RelTol',1e-6,'Stats','on');
x=1;


[tv,Yv{i}]=ode23tb(@(t,y) funsys_2(t,y,x,G,c,M1,mu,m1,q,alpha,sigma,cs),[0 1e20],[a0;e0],options);
Yv{i}(:,1)=Yv{i}(:,1)/149597870700;


end

semilogy(Yv{1}(:,2),Yv{1}(:,1),Yv{2}(:,2),Yv{2}(:,1),Yv{3}(:,2),Yv{3}(:,1),Yv{4}(:,2),Yv{4}(:,1),Yv{5}(:,2),Yv{5}(:,1))
xlabel('Eccentricity')
ylabel('Semi-Major Axis [AU]')
legendCell = cellstr(num2str(q'));
hleg = legend(legendCell,'Location','southeast');
htitle = get(hleg,'Title');
set(htitle,'String','Mass ratio q')
grid on;
axis([0 1 10^(-5) 10])
