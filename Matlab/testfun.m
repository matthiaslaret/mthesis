function testfun

x=linspace(0,1,100);

small = (1-x.^2)./(x + 1./(10000*x)) .* (0.5 - 1./sqrt(1-x.^2));

% x1=linspace(0.2,1,100);
% c5 = 0.9747;
% c4 = -0.7*c5;
% greater = (c4./x1 + c5);

plot([x],[small])
grid on;

end

