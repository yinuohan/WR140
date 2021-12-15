function E = true_to_eccentric(theta, eccentricity)

E = atan2d(sqrt(1 - eccentricity^2) * sind(theta), (eccentricity + cosd(theta)));

end