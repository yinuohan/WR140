% Run spiral model first
path = "C:\\users\\Yinuo\\Desktop\\WR140\\";
table = readtable(path + "Params_new.xlsx");
%average_speeds = table.Windspeed';
average_speeds = table.Bar';
%average_speeds = table.Ellipse';

year = 3.154e+7; % s
AU = 1.496e+8; % km
period = 2896.35/365.25; % yr

%begin_phase = 0;
begin_phase = 0.041;
%begin_phase = -0.041;
% Pos: fitting to bar
% Neg: fitting to ellipse
phases = [0.043, 0.059, 0.077, 0.111, 0.183, 0.592];
total_times = (phases - begin_phase) * period; % yr
segment_times = diff([0 total_times]);
%central_phases = [0 phases(1:end-1)] + segment_times/2;

total_distances = average_speeds .* total_times * year / AU; % AU

segment_distances = diff([0 total_distances]); % AU
segment_speeds = segment_distances ./ segment_times * AU / year;

segment_speeds_diff = diff(segment_speeds);
acceleration = segment_speeds_diff ./ segment_times(2:end);

figure
subplot(2,2,1)
stairs(total_times, segment_speeds, 'o-')
xlabel('Time since dust produced (yr)')
ylabel('Dust speed (km/s)')

subplot(2,2,2)
stairs(total_times(2:end), acceleration, 'o-')
xlabel('Time since dust produced (yr)')
ylabel('Acceleration (km/s / year)')

subplot(2,2,3)
stairs([0 total_times], [0 total_distances], 'o-')
xlabel('Time since dust produced (yr)')
ylabel('Displacement (AU)')

