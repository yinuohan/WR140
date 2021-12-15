function C2 = rotate_y(angle)
    C2 = zeros(3,3);
    C2(:,1) = [cos(angle), 0, -sin(angle)];
    C2(:,2) = [         0, 1,           0];
    C2(:,3) = [sin(angle), 0,  cos(angle)];
end