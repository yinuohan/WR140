% ---------- Data image ----------
% Image files
path = 'C:\\users\\Yinuo\\Desktop\\WR140\\';
epoch1 = [
    "NIRC_ffa_kcont_99-07-29_nodust.fits",
    "NIRC_ffa_pahcs_99-07-29_nodust.fits",
    "NIRC_ffa_feii_99-07-29_nodust.fits"
];
epoch2 = [
    "NIRC_annulus_h_01-06-11.fits",
    "NIRC_annulus_k_01-06-11.fits",
    "NIRC_annulus_ch4_01-06-11.fits"
];
epoch3 = [
    "NIRC_annulus_h_01-07-29.fits",
    "NIRC_annulus_k_01-07-29.fits",
    "NIRC_ffa_oii_01-07-29.fits",
    "NIRC_ffa_kcont_01-07-29.fits",
    "NIRC_ffa_pahcs_01-07-29.fits",
    "NIRC_ffa_feii_01-07-29.fits"
];
epoch4 = [
    "NIRC_ffa_kcont_02-07-23.fits",
    "NIRC_ffa_pahcs_02-07-23.fits"
];
epoch5 = [
    "NIRC_annulus_j_04-05-28_nodust.fits",
    "NIRC_annulus_k_04-05-28_nodust.fits"
];

% Best images % Upper clip
file1 = epoch1(1); % 0.1 no dust
file2 = epoch2(2); % 0.1
file3 = epoch3(2); % 0.01
file4 = epoch4(2); % 0.01
file5 = epoch5(1); % 0.01 no dust

obs_date1 = 2451388.5/365; % 99-07-29
obs_date2 = 2452071.5/365; % 01-06-11
obs_date3 = 2452119.5/365; % 01-07-29
obs_date4 = 2452478.5/365; % 02-07-23
obs_date5 = 2453518.5/365; % 04-05-28

close all

% Pick one epoch
file = file5
obs_date = obs_date5
image = fitsread(path + file);

% Find centre
% [M I] = max(image(:));
% [x y] = ind2sub(size(image),I)

% Cut out
% dim = 128
% lim = floor(dim/2)
% image = image(x-lim:x+lim,y-lim:y+lim);

% Normalise
image = image/max(max(image));
image = min(max(image,0),0.01);

% Log scale
% image = log(image);
% image = min(max(image,-8),-2);

% Renormalise
image = image - min(min(image));
image = image/max(max(image));

% Filter
% blur = imgaussfilt(image,15);
% skeleton = image - blur;
% skeleton = imgaussfilt(skeleton,1);
% skeleton = min(max(skeleton,-0.05),0.15);
% skeleton = skeleton/max(max(skeleton));

% Flip image
%skeleton = image;
skeleton = fliplr(skeleton);

if 0
    figure()
    imagesc(-image)
    set(gca,'YDir','normal')
    axis image
    colormap gray
end

if 1
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
theta_lim = [-120, 160];
turn_off = 0;
n_circ = 2.5

offset = obs_date - n_circ*period - periastron
% This line fixes the above issue with n_circ having to be an integer

% Testing offset
%offset = 7

% Generate Spiral
make_gif = 0;
dim = 128;
pix = 10;

if ~make_gif
    %figure
    [im, theta] = spiral(skeleton,make_gif,dim,pix,windspeed,period,inclination,big_omega,turn_off,eccentricity,omega_lock,little_omega,periastron,cone_angle,offset,n_circ,theta_lim);
    %figure
    %[im, theta] = spiral_idl(skeleton,make_gif,windspeed,period,inclination,big_omega,turn_off,eccentricity,omega_lock,little_omega,periastron,cone_angle,offset,n_circ,theta_lim);
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
    
    save_gif = 0
    
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
