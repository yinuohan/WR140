% ---------- Data image ----------
% Get data
path = "C:\\users\\Yinuo\\Desktop\\WR140\\";

%close all

% Index of image to plot
%indices = [1 2 3 4 9 12 18];

% Order of orbital phase
%indices = [9 12 1 2 18 3 4];

% Grouped by instruments
indices = [1 2 3 4];
indices = [9 12 18];

% Grouped by phase
indices = [9 12 1 2];
indices = [18 3 4];

% Modifications
indices = [12 3];

n_images = length(indices);

figure

for group = 1:n_images

index = indices(group);

% Get plotting parameters
table = readtable(path + "Params.xlsx");

table.Obs_date = table.Obs_date/365.25;

image = fitsread(path + table.Filename(index));

% Find centre
if index == 2
    image = image(250:end, 1:250);
end

[M I] = max(image(:));
[y x] = ind2sub(size(image),I);

x = x + table.Xshift(index);
y = y + table.Yshift(index);

% Cut out
dim = table.Dim(index);
pix = table.Pix(index);

cut_out_these = [1,2,3,4,21];

if sum(index == cut_out_these)
    lim = floor(dim/2);
    image = image(y-lim+1:y+lim,x-lim+1:x+lim);
end

% Normalise
image = image/max(max(image));
image = min(max(image,table.Lowcutoff(index)),table.Highcutoff(index));

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
if table.Flip(index)
    skeleton = fliplr(image);
else
    skeleton = image;
end


% ---------- Model ----------
% Fixed
omega_lock = false;
period = 2896.35/365.25;
eccentricity = 0.8964;
little_omega = 46.8;
periastron = 2446155.3/365.25;
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
windspeed = 2350 * 0.210805 / distance; % km/s to mas/year
% windspeed = 2460 -130 +130

cone_angle = 38*2;
% cone angle = (42 -5 +5) * 2

theta_lim = [-138 135]; 
% Low = -(135 -5 +10), High = 135 -5 +15


% On and off
turn_off = 0;
n_circ = 1;

obs_date = table.Obs_date(index);
offset = obs_date - n_circ*period - periastron;
% The programs implies: {n_circ} periods ago, the position angle of the
% binary corresponds to an offset from periastron of {offset} years


% Generate Spiral
make_gif = 0;
save_gif = 0;

phase = (table.Obs_date(index) - periastron)/period;
phase = round(phase - floor(phase), 2);

Title = string(table.Filter(index)) + " \phi = " + string(phase) + " (" + string(table.Time(index)) + ")";

if ~make_gif
    [im, comment] = spiral(skeleton,make_gif,Title,dim,pix,windspeed,period,inclination,big_omega,turn_off,eccentricity,omega_lock,little_omega,periastron,cone_angle,offset,n_circ,theta_lim);
    
    scale = dim/1e3;
    dim = pix/2;
    x = linspace(-dim*scale, dim*scale, dim*2);
    
    % ---------- PLOT DATA ONLY ----------
    ax(3*(group-1)+1) = subplot(n_images, 3, 3*(group-1)+1);
    imagesc(-skeleton)
    if group == length(indices)
        xlabel("Relative RA ('')")
    end
    ylabel("Relative Dec ('')")
    set(gca,'YDir','normal')
    axis image
    %colormap gray
    title(Title)
    colormap(ax(3*(group-1)+1), gray)
    
    % ---------- PLOT OVERLAY ----------
    % DATA
    ax(3*(group-1)+2) = subplot(n_images, 3, 3*(group-1)+2);
    imagesc(-x,x,skeleton)
    
    % MODEL OUTLINE
    hold on
    high_pass_model = imgaussfilt(im - imgaussfilt(im,2),1);
    model_threshold = 0.0005;
    high_pass_model(high_pass_model > model_threshold) = 1;
    high_pass_model(high_pass_model <= model_threshold) = 0;
    high_pass_model = imgaussfilt(high_pass_model,1);
    high_pass_model = high_pass_model/max(max(high_pass_model));
    
    h2 = imagesc(-x,x,high_pass_model);
    if group == length(indices)
        xlabel("Relative RA ('')")
    end
    %ylabel("Relative Dec ('')")
    axis image
    %colormap gray
    %colormap parula
    set(gca,'YDir','normal')
    set(gca,'XDir','reverse')
    alpha = (1 - (high_pass_model<0.85)) * 0.8;
    set(h2, 'AlphaData', alpha)
    %title(Title)
    colormap(ax(3*(group-1)+2), parula)
    
    % ---------- PLOT MODEL ONLY ----------
    ax(3*(group-1)+3) = subplot(n_images, 3, 3*(group-1)+3);
    %imagesc(x,x,-min(im,0.1));
    imagesc(-min(im,0.1));
    if group == length(indices)
        xlabel("Relative RA ('')")
    end
    %ylabel("Relative Dec ('')")
    %colormap gray
    axis image
    set(gca,'YDir','normal')
    title(comment)
    colormap(ax(3*(group-1)+3), gray)
    
    %colormap(ax(1), gray)
    %colormap(ax(2), parula)
    %colormap(ax(3), gray)
end

end


% ---------- Animation ----------
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
