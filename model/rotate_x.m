function C2 = rotate_x(angle)
    C2 = zeros(3,3);
    C2(:,1) = [1, 0, 0];
    C2(:,2) = [0,  cos(angle), sin(angle)];
    C2(:,3) = [0, -sin(angle), cos(angle)];
end