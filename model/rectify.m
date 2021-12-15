function theta2 = rectify(theta)
    x = theta/(2*pi);
    theta2 = 2*pi*(x + floor(0.5-x));
end