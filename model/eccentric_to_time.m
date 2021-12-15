function time = eccentric_to_time(E, period, eccentricity)

time = period/2/pi * (E/180*pi - eccentricity * sind(E));

end