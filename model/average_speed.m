% Run spiral model first
path = "C:\\users\\Yinuo\\Desktop\\WR140\\";
table = readtable(path + "Params_new.xlsx");
average_speeds = table.Windspeed';

begin_phase = 0;
phases = [0.043, 0.059, 0.077, 0.111, 0.183, 0.592] - begin_phase;
total_times = phases;
total_distances = average_speeds .* total_times;

segment_times = diff([begin_phase total_times]);
segment_distances = diff([0 total_distances]);
segment_speeds = segment_distances ./ segment_times;

segment_speeds_diff = diff([0 segment_speeds]);
acceleration = segment_speeds_diff ./ segment_times;

figure
subplot(2,2,1)
stairs(total_times, segment_speeds, 'o-')
xlabel('Orbital phase')
ylabel('Dust speed (km/s)')

subplot(2,2,2)
stairs(total_times, acceleration, 'o-')
xlabel('Orbital phase')
ylabel('Acceleration (km/s / year)')

subplot(2,2,3)
stairs([0 total_times], [0 total_distances], 'o-')
xlabel('Orbital phase')
ylabel('Displacement (km)')

% figure
% subplot(2,1,1)
% stairs(total_distances, segment_speeds, 'o-')
% xlabel('Orbital phase')
% ylabel('Dust speed (km/s)')
% 
% subplot(2,1,2)
% stairs(total_distances, acceleration, 'o-')
% xlabel('Orbital phase')
% ylabel('Acceleration (km/s / year)')

