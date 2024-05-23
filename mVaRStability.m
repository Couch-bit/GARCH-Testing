x=linspace(-2.5, 2.5, 1000);
y=linspace(0, 12.5, 1000);
[X,Y] = meshgrid(x, y);
ineq = (abs(X) <= 6*(sqrt(2) - 1)) ...
& (27*Y.^2 - (216 + 66*X.^2).*Y + 40*X.^4 + 336*X.^2 <= 0);
ineq_fix = double(ineq);
ineq_fix(ineq_fix==0) = NaN;
h = pcolor(X,Y,double(ineq_fix));
h.EdgeColor = 'none';
title("Dziedzina Dobrej Określoności dla mVaR")
xlabel("Skośność");
ylabel("Eksces"); 
colormap winter;
grid on;