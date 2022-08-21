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
%speed = table.Windspeed(index);
speed = 2400; % km/s
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
n_circ = 2;

% The program implies: {n_circ} periods ago, the position angle of the
% binary corresponds to an offset from periastron of {offset} years


% Generate Spiral
make_gif = 1;
save_gif = 1;

if make_gif
    h = figure;
    axis tight manual % this ensures that getframe() returns a consistent size
    filename = '/users/Yinuo/Desktop/WR140.gif';
    
    scale = 10;
    pix = scale;
    dim = 512;
    x = linspace(-dim/2*scale/1000, dim/2*scale/1000, dim);
    
    offset_initial = 0
    increment = 0.1
    offset_final = floor(period*10)/10
    
    for offset = offset_initial:increment:offset_final
        model = spiral(0,make_gif,'',dim,pix,windspeed,period,inclination,big_omega,turn_off,eccentricity,omega_lock,little_omega,periastron,cone_angle,offset,n_circ,theta_lim);
        
        hold off
        limit = 500;
        imagesc(x,x,min(imgaussfilt(model,1),limit));
        title("Offset: " + string(offset) + " yr")
        xlabel("Relative RA ('')")
        ylabel("Relative Dec ('')")
        %colormap gray
        colormap hot
        axis image
        set(gca,'YDir','normal')
        set(gcf,'Position',[20 20 580 580])

        hold on
        plot([-2,-1], [-2.3,-2.3], '-w', 'LineWidth',1)
        text(-1.72,-2.15,"1 arcsec",'Color','w')
        text(1.5,-2.2,"Offset (yr): " + string(offset),'Color','w')
        pause(0.01)

        if save_gif
            if offset == offset_initial
                gif(filename, 'nodither')
            else
                gif('nodither', 'frame', gca)
            end
        end

    end
end
