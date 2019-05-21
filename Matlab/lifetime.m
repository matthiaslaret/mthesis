function [lifetime]=lifetime(x,dec,M,mu,a0,e0,R,H,alpha,Md,sigma,cs)
%eccentricity evolution for different initial eccentricities with fixed
%initial semi-major axis


%constants
G = 6.67408e-11;
c = 299792458;
% G = 1;
% c = 1;



close all;

cc=jet(30);
i = 1;


options = odeset('RelTol',1e-6,'Stats','on');
if dec == 1
options = odeset('RelTol',1e-5,'events',@events);
end
if dec == 2
options = odeset('RelTol',1e-5,'events',@events);
end
if dec == 3
options = odeset('RelTol',1e-5,'events',@events);
end
[tv,Yv]=ode23tb(@(t,y) funsys_2(t,y,x,G,c,M,mu,H,R,alpha,Md,sigma,cs),[0 1e20],[a0;e0],options);
txt = ['Initial e_0  ',num2str(e0)];

if x==0
tv = tv*10^12;
end


plot(tv/31536000,Yv(:,1),'-','color',cc(i,:),'DisplayName',txt);
plot(tv/31536000,Yv(:,2),'-','color',cc(i,:),'DisplayName',txt);
lifetime = tv(end)/31536000;
if dec==1
    lifetime = tv(end)/31536000;
end
if dec==2
    lifetime = Yv(end,2);
end
if dec==3
    lifetime = Yv(end,1);
end
xlabel({'Time','years'})
ylabel('Eccentricity')


hold off
legend show
grid
 

function [value,isterminal,direction] = events(t,Y)
  value = (G^3/c^5 * 64/5 * (mu*M^2)/Y(1)^3 * 1/(1-Y(2)^2)^(7/2) *( 1 + 73/24 * Y(2)^2 + 37/96 * Y(2)^4)) - 24*pi*alpha*Y(1)*sigma*cs^2 / (mu*sqrt(G*M/(Y(1)^3))); 
  isterminal = 1;
  direction = 0;
end




%semi-major axis evolution for different initial semimajor axis

% 
% cc=jet(30);
% i = 1;
% 
% hold on
% for a0 = 1:1:20
%       
% e0 = 0.8;
% 
% [tv,Yv]=ode45('funsys',[0 200],[a0;e0]);
% txt = ['Initial a_0  ',num2str(a0)];
% plot(tv,Yv(:,1),'-','color',cc(i,:),'DisplayName',txt);
% xlabel({'Time','units?'})
% ylabel('Semi-Major Axis')
% 
% i = i+1;
% 
% end
% 
% hold off
% legend show
% grid


% % test with gas
% figure;
% a0 = 10;
% e0 = 0.6;
% [tv,Yv]=ode45('funsys',[0 200],[a0;e0]);
% [tv_gas,Yv_gas]=ode45('funsys_gas',[0 200],[a0;e0]);
% subplot(2,1,1);
% plot(tv,Yv(:,1),'.-',tv_gas,Yv_gas(:,1),'--');
% ylabel('Semi-Major Axis');
% legend('GW','GW + gas','location','northeast','Orientation','horizontal');
% 
% 
% subplot(2,1,2);
% plot(tv,Yv(:,2),'.-',tv_gas,Yv_gas(:,2),'--');
% ylabel('Eccentricity');
% 
% suptitle(['a_0 = ',num2str(a0),'; e_0 = ',num2str(e0)]);

% %a(e)
% 
% e=0:0.01:1;
% a = e.^(12/19) .* (1./(1-e.^2)) .* (1 + 121/304 * e.^2).^(870/2299);
% plot(e,a)
end