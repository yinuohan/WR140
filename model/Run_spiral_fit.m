% ---------- Data image ----------
close all
path = "/Users/yinuo/Desktop/WR140/WR140/";
outdir = "/Users/yinuo/Desktop/WR140 Nature Rev2/Figures/Models Ellipse Rev2/";

% Get data
%figure('Renderer', 'painters', 'Position', [100 100 1200 900])
%feature = "bar";
feature = "ellipse";
index = 13;
suffix = "";
%suffix = "_noacc";

% Get plotting parameters and image
table = readtable("/Users/yinuo/Desktop/WR140/WR140 nature rev 1/Figures nature rev 1/Params_rev1.xlsx");
table.Obs_date = table.Obs_date/365.25;
image = fitsread(path + "NIRC2_stacked/" + table.Filename(index));

% Normalise
image = image/max(max(image));
image = min(max(image,table.Lowcutoff(index)),table.Highcutoff(index));
image = image/max(max(image));


% ---------- Model ----------
% Fixed
omega_lock = false;
period = 2896.35/365.25;
eccentricity = 0.8964;
little_omega = 46.8;
periastron = 2446155.3/365.25;
big_omega = 353.6;
%inclination = 119.6;
inclination = 180;
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

% Adjustable
%speed = table.Windspeed(index);
if feature == "ellipse"
    speed = table.Ellipse(index);
elseif feature == "bar"
    speed = table.Bar(index);
else
    disp("Unrecognised feature")
end
windspeed = speed * 0.21095 / distance; % km/s to mas/year
% windspeed = 2460 -130 +130
% 2400

cone_angle = 40*2;
% cone angle = (42 -5 +5) * 2
% 40
% 38

theta_lim = [-135 135]; 
% Low = -(135 -5 +10), High = 135 -5 +15

% On and off
turn_off = 0;
n_circ = 1;

obs_date = table.Obs_date(index);
offset = obs_date - n_circ*period - periastron;
% The program implies: {n_circ} periods ago, the position angle of the
% binary corresponds to an offset from periastron of {offset} years


% Generate Spiral
make_gif = 0;
save_gif = 0;

phase = (table.Obs_date(index) - periastron)/period;
phase = round(phase - floor(phase), 3);

dim = 512;
datadim = table.Dim(index);
scale = dim/datadim;
datapix = table.Pix(index);
pix = datapix/scale;

datax = linspace(-dim*pix/2/1e3, dim*pix/2/1e3, datadim);
x = linspace(-dim*pix/2/1e3, dim*pix/2/1e3, dim);

Title = " \phi = " + string(phase);
[model, comment] = spiral(image,make_gif,Title,dim,pix,windspeed,period,inclination,big_omega,turn_off,eccentricity,omega_lock,little_omega,periastron,cone_angle,offset,n_circ,theta_lim);
im = model;

% ---------- PLOT OVERLAY ----------
imagesc(-datax,datax,image)

% Blur to make same as Python version
%model = imgaussfilt(model,2);

orig_model = model;
model_threshold = 0.1 * max(max(model));
model = min(max(model, 0), model_threshold);
model = imgaussfilt(model,2);
model = model - imgaussfilt(model,2);
if index == 13 || index == 14
    model_threshold2 = 0.05 * max(max(model));
else
    model_threshold2 = 0.01 * max(max(model));
end
model = (model >= model_threshold2);
%alpha = (model > model_threshold2) * 0.8;
alpha = (model > 0.1) * 0.7;
model = 10 * model / max(max(model));

hold on
h2 = imagesc(-x,x,model);
axis image
set(gca,'YDir','normal')
set(gca,'XDir','reverse')
set(h2, 'AlphaData', alpha)
title(Title)
caxis([0 1])
colormap parula

%fitswrite(im, outdir + "Rev1_" + string(phase) + suffix + ".fits")
%fitswrite(im, outdir + "Rev2_" + string(phase) + suffix + ".fits")
fitswrite(im, outdir + "Rev2VISorbital_azimuthal_" + string(phase) + suffix + ".fits")

%imagesc(orig_model)
%set(gca,'YDir','normal')
%set(gca,'XDir','normal')
