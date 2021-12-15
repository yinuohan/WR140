function eta = eta(theta)

theta = theta/180*pi;

eta = (tan(theta) - theta) / (tan(theta) - theta + pi);

end