

x = [10 20 30 40 50 60 70 80];
y = [25 70 380 550 610 1220 830 1450];
p=linregr(x,y);
p=[p(1) p(2)];
yestimation = polyval(p, x)
figure; 
plot(y,'r');
hold on;
plot(yestimation,'b');
grid on;