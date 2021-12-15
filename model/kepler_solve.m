function theta = kepler_solve(t, P, ecc)
    maxj = 50; % Max number of iterations
    tol = 1e-8; % Convergence tolerance
    
    M = 2*pi/P*t;
    E = zeros(1,length(t));
    %tj = 0;
    
    for i = 1:length(t)
        E0 = M(i);
        %E1s = zeros(1,length(t))
        
        % Newton's formula to solve for eccentric anomoly
        for j = 1:maxj
            E1 = E0 - (E0 - ecc*sin(E0)-M(i))/(1 - ecc*cos(E0));
           %E1s(j) = E1;
            if abs(E1 - E0) < tol
                break
            end
            E0 = E1;
        end
        
        if j == maxj
            disp('Did not converge')
        end
        
        E(i) = E1;
        %tj = tj+j;
    end
    
    %plot(E1s)
    
    % Compute 2-dimensional spiral angles & radii
    theta = 2*atan(sqrt((1+ecc)/(1-ecc))*tan(E/2));
end