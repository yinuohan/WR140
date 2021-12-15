function radius = true_to_radius(theta, a, e)

radius = a * (1 - e^2) ./ (1 + e * cosd(theta));

end
