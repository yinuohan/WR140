% ---------- Data image ----------
% Image files
path = 'C:\\users\\Yinuo\\Desktop\\WR140\\';

%close all

files = [
    "n0176_2017_7_31.fits",
    "n0126_drp_2017_11_5.fits",
    "WR140_NIRC2_Lp.fits",
    "WR140_NIRC2_Ms.fits"
    ];
obs_dates = [
    2457965,
    2458062,
    2453664,
    2453664
    ];
lowcutoff = [
    0.1,
    0,
    0,
    0
    ];
highcutoff = [
    1,
    0.8,
    0.01,
    0.01
    ];
yshift = [
    5,
    0,
    0,
    0
    ];
xshift = [
    2,
    10,
    0,
    0
    ];
dims = [
    128,
    128,
    420,
    420
    ];
pixs = [
    9.952,
    9.952,
    9.952,
    9.952
    ];

index = 3;

file = files(index);
obs_date = obs_dates(index)/365;

image = fitsread(path + file);

% Find centre
if index == 2
    image = image(250:end, 1:250);
end

[M I] = max(image(:));
[y x] = ind2sub(size(image),I);

x = x + xshift(index);
y = y + yshift(index);

% Cut out
dim = dims(index);
pix = pixs(index);

if index == 1 || index == 2 || index == 3 || index == 4
    lim = floor(dim/2);
    image = image(y-lim+1:y+lim,x-lim+1:x+lim);
end

% Normalise
image = image/max(max(image));
image = min(max(image,lowcutoff(index)),highcutoff(index));

% Log scale
% image = log(image);
% image = min(max(image,-8),-2);

% Renormalise
image = image/max(max(image));

% Filter
% blur = imgaussfilt(image,15);
% skeleton = image - blur;
% skeleton = imgaussfilt(skeleton,1);
% skeleton = min(max(skeleton,-0.05),0.15);
% skeleton = skeleton/max(max(skeleton));

% Flip image
skeleton = image;

if 0
    figure()
    imagesc(-image)
    set(gca,'YDir','normal')
    axis image
    colormap gray
end

if 0
    figure()
    imagesc(-skeleton)
    set(gca,'YDir','normal')
    axis image
    colormap gray
end

% ---------- Model ----------
% Fixed
omega_lock = false;
period = 2896.35/365;
eccentricity = 0.8964;
little_omega = 46.8;
periastron = 2446155.3/365;
big_omega = 353.6;
inclination = 119.6;
distance = 1.67;

% Masses
mwc = 14.9; % ± 0.5 Msum for the WC7
mo =  35.9; % ± 1.3 Msun for the O5
M = mwc + mo; % Enclosed mass
a = (M*period^2)^(1/3); % semimajor axis in AU

% Constants
Msun = 1.98847e30; % kg
G = 6.67408e-11; % SI units
AU = 1.495978707e11; % m

% Vis-Viva
d = a*(1-eccentricity); % AU
v = sqrt( G*M*Msun * (2/d/AU - 1/a/AU) ) *1e-3; % km/s

% big_omega = 0
% little_omega = 0
% inclination = 0
% eccentricity = 0

% Adjustable
windspeed = 2460 * 0.210805 / distance; % km/s to mas/year
cone_angle = 43*2;
%offset = -1
%obs_date = 2452071.50000/365 % JD of 01-06-11
% offset = obs_date - periastron
% In order for offset to correspond to current position angle, 
% n_circ must be integer
% The programs implies: {n_circ} periods ago, the position angle of the
% binary corresponds to an offset from periastron of {offset} years

% On and off
theta_lim = [-140 140]; % [-130 160] [-140 140]
turn_off = 0;
n_circ = 1;

offset = obs_date - n_circ*period - periastron;
% This line fixes the above issue with n_circ having to be an integer

% Testing offset
%offset = 0

% Generate Spiral
make_gif = 0;
save_gif = 0;

if ~make_gif
    %figure
    [im, theta] = spiral(skeleton,make_gif,dim,pix,windspeed,period,inclination,big_omega,turn_off,eccentricity,omega_lock,little_omega,periastron,cone_angle,offset,n_circ,theta_lim);
end

if ~make_gif && 0
    figure
    plot(theta)
end

if make_gif
    h = figure;
    axis tight manual % this ensures that getframe() returns a consistent size
    filename = 'C:/users/Yinuo/Desktop/WR140.gif';
    
    scale = 45/1e3;
    dim = 360/2;
    x = linspace(-dim*scale, dim*scale, dim*2);
    
    offset_initial = 0
    increment = 0.2
    offset_final = ceil(period)
    
    for offset = offset_initial:increment:offset_final
        model = spiral(skeleton,make_gif,dim,pix,windspeed,period,inclination,big_omega,turn_off,eccentricity,omega_lock,little_omega,periastron,cone_angle,offset,n_circ,theta_lim);
        
        hold off
        limit = 50;
        imagesc(x,x,-min(model,limit));
        title("Offset: " + string(offset) + " yrs")
        xlabel("Relative RA ('')")
        ylabel("Relative Dec ('')")
        colormap gray
        axis image
        set(gca,'YDir','normal')
        set(gcf,'Position',[20 20 580 580])

        hold on
        plot([-6,-4], [-7.3,-7.3], '-k', 'LineWidth',1)

        text(-6,-6.8,"2 arcsec")
        text(-6,-6.1,"Offset (yr): " + string(offset))
        pause(0.01)

        if save_gif
            if offset == offset_initial
                gif(filename, 'nodither')
            else
                gif('nodither', 'frame',gca)
            end
        end

    end
end
