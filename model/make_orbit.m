% Input
a = 100
e = 0.9

omega = 0
inc = 0
Ohm = 0

% Parameters
c = e * a;
b = sqrt(a^2 - c^2);

%ellipse = ellipse * rotate_z(-omega) * rotate_y(-inc) * rotate_z(-Ohm);

% Plot ellipse
z = 0;

theta = 0:1:360;
x = a*cosd(theta);
y = b*sind(theta);

figure
hold on
plot(x, y)
scatter(0, 0, '+')
scatter(c, 0, '+')
axis equal

% Add orbit
period = 100;
t = 0:1:2*period;
angle = kepler_solve(t, period, e)/pi*180;
x = a*cosd(angle);
y = b*sind(angle);

for i = 1:length(t)
    scatter(x(i), y(i), 'r.')
    drawnow
end

figure
hold on
plot(angle)
plot(diff(angle))