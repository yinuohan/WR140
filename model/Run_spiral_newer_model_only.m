% ---------- Data image ----------
% Get data
path = "C:\\users\\Yinuo\\Desktop\\WR140\\";
outdir = "C:\\users\\Yinuo\\Desktop\\WR140\\NIRC2_stacked\\";

%close all

% Times to generate
obs_dates = [2452071.5, 2452119.5, 2457965, 2458062, 2452478.5, 2453664]/365.25;
outnames = ["01jun", "01jul", "17jul", "17nov", "02jul", "05oct"];
phases = [0.0426, 0.0592, 0.0774, 0.1109, 0.1832, 0.5925];
dims = [128, 128, 128, 128, 128, 512];
pixs = [10, 10, 9.52, 9.52, 10, 9.52];

% Params for high
dims = [512, 512, 512, 512, 512, 512];
pixs = [10/4, 10/4, 9.52/4, 9.52/4, 10/4, 9.52];

for group = 1:length(obs_dates)
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

    % Adjustable
    table = readtable(path + "Params_new.xlsx");
    speed = table.Windspeed(group);
    windspeed = speed * 0.210805 / distance; % km/s to mas/year
    % windspeed = 2460 -130 +130

    cone_angle = 38*2;
    % cone angle = (42 -5 +5) * 2

    theta_lim = [-135 135]; 
    % Low = -(135 -5 +10), High = 135 -5 +15

    % On and off
    turn_off = 0;
    n_circ = 1;

    obs_date = obs_dates(group);
    offset = obs_date - n_circ*period - periastron;
    % The programs implies: {n_circ} periods ago, the position angle of the
    % binary corresponds to an offset from periastron of {offset} years

    % Calculate phase
    phase = (obs_date - periastron)/period;
    phase = round(phase - floor(phase), 3);

    % Generate Spiral
    dim = dims(group);
    pix = pixs(group);
    make_gif = 0;
    Title = "";
    [im, comment] = spiral([],make_gif,Title,dim,pix,windspeed,period,inclination,big_omega,turn_off,eccentricity,omega_lock,little_omega,periastron,cone_angle,offset,n_circ,theta_lim);

    % Write to fits
    disp(phase)
    fitswrite(im, outdir + "MODELhighmodwind_" + outnames(group) + ".fits")
end
