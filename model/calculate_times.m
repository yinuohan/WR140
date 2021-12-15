offset = 20
range = 280

neg = -range/2+offset
pos = range/2+offset

negt = eccentric_to_time(true_to_eccentric(neg, 0.9), 8, 0.9)
post = eccentric_to_time(true_to_eccentric(pos, 0.9), 8, 0.9)
time = post - negt