% Линейная интерполяция радиуса резонатора (письмо от 22.09.23)
Rr = load('d:\Alex\Documents\Work\novozhilova\22-09-23 Продольная структура и профиль рез-ра\RR2812.dat');

plot(Rr(:,1),Rr(:,2))
hold on

L = 10;
x = 0:0.0001:L;
y = zeros(length(x),1);

[b, c, d] = spline(length(Rr(:,1)), Rr(:,1), Rr(:,2));
s = zeros(length(x),1);

x = x/L*(Rr(end,1) - Rr(1,1)) + Rr(1,1);

for i=1:length(x)
   y(i) = uval(x(i),Rr(:,1),Rr(:,2));
   s(i) = seval(x(i), length(Rr(:,1)), Rr(:,1), Rr(:,2), b, c, d);
end


plot(x,s,x,y)
