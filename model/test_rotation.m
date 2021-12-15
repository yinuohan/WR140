x = [1; 0; 0];
y = [0; 1; 0];
z = [0; 0; 1];

rx = rotate_x(pi/2);
ry = rotate_y(pi/2);
rz = rotate_z(pi/2);

a = [0 0 0;
     0 0 0;
     1 0 0];

imagesc(a)
set(gca,'YDir','normal')

