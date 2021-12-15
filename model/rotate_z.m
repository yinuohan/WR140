function C2 = rotate_z(angle)
    C2 = zeros(3,3);
    C2(:,1) = [ cos(angle), sin(angle), 0];
    C2(:,2) = [-sin(angle), cos(angle), 0];
    C2(:,3) = [0, 0, 1];
end