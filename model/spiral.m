function [im, comment] = spiral(skeleton, gif, Title, dim, pix, windspeed,period,inclination,big_omega,turn_off,eccentricity,omega_lock,little_omega,periastron,cone_angle,offset,n_circ,theta_lim)

% Plotting parameters
n_p = 500;                      % Number of points per circle
n_c = 5000;                     % Number of circles per period

im_siz = dim;        %360       % Image size (pixels) (square)
im_res = pix;                   % Image resolution (mas/pix)

% Physics parameters
w     = windspeed / 365.25;        % windspeed [mas/year]->[mas/day]
P     = period * 365.25;           % period [year]->[day]
omega = deg2rad(little_omega);  % omega
inc   = deg2rad(inclination);   % inclination  
Ohm   = deg2rad(big_omega);     % Omega 
ecc   = eccentricity;           % Eccentricity
pa    = periastron * 365.25;       % Periastron date (Julian days)
cone  = deg2rad(cone_angle);    % cone angle (deg)
offst = offset * 365.25;           % time offset [year]->[day]
%offst = t_obs - pa;            % time offset
rnuc = turn_off;
lim = theta_lim/180*pi;         % limit for theta to produce dust

% Don't use in general
if omega_lock
    adjust = kepler_solve(pa-245, P, ecc);
    omega = rectify(omega + adjust);
end

% Time vector
t = (0:n_circ*n_c)/n_c*P + offst; 

% Angles from Keplers laws as a function of time
theta = kepler_solve(t, P, ecc); 

% Radius of each dust circle
r2 = w * (n_circ*P - (0:n_circ*n_c)/n_c*P);
% First ring in time has biggest radius

% 3D spiral plume
% NOTE: Spiral is always generated to align with first point
% on the x-axis (theta = 0). The result is then rotated by 
% the specified angles using direction cosine matrices 
% (note: order of rotations is important & must be preserved
% - omega, inc, Omega)

% Generate the coordinates of the initial unit circle
% which is copied around the spiral
chi = (0:n_p-1)/n_p*pi*2; % Angle
ones2 = ones(1,length(chi)); % Just ones

% The circle is parallel to the y-z plane
circ = [(cos(cone/2)*ones2); % x - becomes North
        (sin(cone/2)*cos(chi)); % y - becomes East
        (sin(cone/2)*sin(chi))]; % z - becomes (-) line of sight

% With anisotropic winds, can try and ellipse
% elongation_factor = 1;
% circ = [(cos(cone/2)*ones2); % x
%         (sin(cone/2)*cos(chi)); % y
%         (sin(cone/2)*sin(chi)) * elongation_factor]; % z

% Initialise full array to store coordinates of all points
circfull = zeros(3, (n_p+1)*length(theta));
gen = 1:n_p; % indices of one circle

% Calculate coordinates of each circle on the spiral
for j = 1:length(theta)
    if r2(j) >= rnuc 
        if rectify(theta(j)) >= lim(1)
            if rectify(theta(j)) <= lim(2)
                circj = rotate_z(theta(j)) * circ * r2(j);
                circfull(:, (j-1)*n_p+gen) = circj;
            end
        end
    end
end

% TEST orbit
% Input
% im_res = 1
% e = ecc
% a = 100
% c = e * a;
% b = sqrt(a^2 - c^2);
% 
% circfull = zeros(3, (n_p+1)*length(theta));
% for theta = 1:1:360
%     circfull(1,floor(theta)) = a*cosd(theta);
%     circfull(2,floor(theta)) = b*sind(theta);
% end
% 
% circfull(1,361) = c;
% circfull(2,361) = 0;
% 
% for i = 362:1:362+a-1
%     circfull(1,floor(i)) = i-362;
%     circfull(2,floor(i)) = 0;
% end

% Rotate points by specifed rotations -------------------
circfull = rotate_z(Ohm) * (rotate_x(inc) * (rotate_z(omega) * circfull));
%circfull = circfull' * rotate_z(omega) * rotate_x(inc) * rotate_z(Ohm) ;
%circfull = circfull';

% Variable density across spiral
use_density = 1; % Varies across angle within each circle
use_d2 = 1; % Varies across spiral so circles differ

comment = "";

if use_density
    % Goes sinusoidally between 1 and 3, reaching 3 at {densest_angle}
    %densest_angle = 230;
    %density = cos(chi - densest_angle/180*pi) + 1; 
    
    % Gaussian
    densest_angle = 180;
    % densest angle = 230 -30 +20
    % densest_angle = 180 fixed, used in paper
    sd = 80;
    % sd = 60 -10 +30;
    % sd = 80 -20 +20;
    dist_chi = min([abs(chi/pi*180-densest_angle); abs(chi/pi*180+360-densest_angle); abs(chi/pi*180-360-densest_angle)]);
    density = exp( - ( (dist_chi/180*pi) / (sd/180*pi) ).^2 );
    
    comment = comment + " " + string(densest_angle) + "/" + string(sd);
end

if use_d2
    % Turn off completely at some stage
    %d2 = (theta < -120/180*pi) + (theta > 120/180*pi);
    
    % Inverse Gaussian
    %least_dense = 0; % perisatron
    %d2 = 1 - exp(-((theta - least_dense/180*pi) / (50/180*pi)).^2);
    
    % Double Gaussian
    dense1 = -135;
    dense2 = 135;
    sd = 40;
    % sd = 40 -10 +30
    
    dist_theta1 = min([abs(theta/pi*180-dense1); abs(theta/pi*180+360-dense1); abs(theta/pi*180-360-dense1)]);
    dist_theta2 = min([abs(theta/pi*180-dense2); abs(theta/pi*180+360-dense2); abs(theta/pi*180-360-dense2)]);
    
    d2 = exp( - ( (dist_theta1/180*pi) / (sd/180*pi) ).^2 ) + exp( - ( (dist_theta2/180*pi) / (sd/180*pi) ).^2 );

    comment = comment + " " + string(dense1) + "/" + string(dense2) + "/" + string(sd);
    
    % Plotting
%     distances = true_to_radius(theta/pi*180, 1, eccentricity);
%     phases = mod(t/365.25.25, period)/period;
%     
%     d2_plot = d2;
%     d2_plot(theta/pi*180 < -135) = 0;
%     d2_plot(theta/pi*180 > 135) = 0;
%     
%     figure
%     hold on
%     %plot(theta)
%     plot([phases 1+phases], [distances distances], '.')
%     plot([phases 1+phases], [d2_plot d2_plot], '.')
%     
%     figure
%     hold off
%     plot(distances, d2_plot)
%     figure
end

% Generate image
n_total = n_p*length(theta);
im = zeros(im_siz, im_siz);

% TEST line
% circfull = zeros(3, (n_p+1)*length(theta));
% for i = 1:100
%     circfull(1,i) = i;
%     circfull(2,i) = i/2;
% end

% Project 3D spiral to pixel values by looping over all points
for i = 1:n_total
    
    % View along z axis
    imy = fix(circfull(1, i)/im_res + im_siz/2);
    imx = im_siz - fix(circfull(2, i)/im_res + im_siz/2);
    
    %imy = fix(circfull(1, i)/im_res + im_siz/2);
    %imx = fix(circfull(2, i)/im_res + im_siz/2);
    
    % Image y = up = North = coordinate x
    % Image -x = left = East = coordinate y
    
    % Add density
    if (imx > 0) && (imx <= im_siz) && (imy > 0) && (imy <= im_siz)
        if use_density && ~use_d2
            dens = density(mod(i-1,n_p)+1);
            im(imy, imx) = im(imy, imx) + dens;
        elseif use_d2 && ~use_density
            dens2 = d2(floor((i-1)/n_p)+1);
            im(imy, imx) = im(imy, imx) + dens2;
        elseif use_d2 && use_density
            dens3 = density(mod(i-1,n_p)+1) * d2(floor((i-1)/n_p)+1);
            im(imy, imx) = im(imy, imx) + dens3;
        else
            im(imy, imx) = im(imy, imx) + 1;
        end
    % Coordinate 1 = image y
    % Coordinate 2 = image x
    end
end

% Get rid of centre
%im(im_siz/2, im_siz/2) = 0;
im(im_siz/2, im_siz/2) = ( sum(sum( im(im_siz/2-1:im_siz/2+1, im_siz/2-1:im_siz/2+1) )) - im(im_siz/2, im_siz/2))/8;

% Normalise
if max(max(im)) > 0 && ~gif && 0
    im = im/max(max(im));
end


% ---------- PLOTTING ----------
scale = im_res/1e3;
dim = im_siz/2;
x = linspace(-dim*scale, dim*scale, dim*2);

if ~gif && 0
    
    % DATA
    figure
    imagesc(x,x,skeleton)
    axis image
    xlabel("Relative RA ('')")
    ylabel("Relative Dec ('')")
    set(gca,'YDir','normal')
    
    % MODEL OUTLINE
    hold on
    high_pass_model = imgaussfilt(im - imgaussfilt(im,2),1);
    model_threshold = 0.0005;
    high_pass_model(high_pass_model > model_threshold) = 1;
    high_pass_model(high_pass_model <= model_threshold) = 0;
    high_pass_model = imgaussfilt(high_pass_model,1);
    high_pass_model = high_pass_model/max(max(high_pass_model));
    
    h2 = imagesc(x,x,high_pass_model);
    xlabel("Relative RA ('')")
    ylabel("Relative Dec ('')")
    axis image
    %colormap gray
    set(gca,'YDir','normal')
    alpha = (1 - (high_pass_model<0.85)) * 0.8;
    set(h2, 'AlphaData', alpha)
    title(Title)
end

% MODEL
if ~gif && 0
    figure
    %h = imagesc(x,x,-min(im,0.1));
    imagesc(-min(im,0.1));
    xlabel("Relative RA ('')")
    ylabel("Relative Dec ('')")
    colormap gray
    axis image
    set(gca,'YDir','normal')
    title(Title)
end


