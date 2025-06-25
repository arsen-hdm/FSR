% This script is used to extract the track references, the medium trajectory
% and the optimal one for racing from
% .csv files

% Of course there are different ways to obtain track's data and racelines,
% like for example for the medium one we also use triangulation in this
% project, but we have to mention the site from which we took the data:
% https://github.com/TUMFTM/racetrack-database/tree/master

clc;clear all; close all;

path="circuit_IMS\";

% In the "track_name".csv file there's already the medium trajectory
% between the inner and the outher borders, so here we're extracting them to
% save firstly the non-optimal trajectory

data = readtable(path+'track.csv');
x = data.x_m;
y = data.y_m;

xp = x';
yp = y';
Pp = [xp, yp];

%save(path+'medium_trajectory.mat', 'Pp', 'xp', 'yp');

% Here instead we're retrieving informations regarding the borders of the
% track, so to have a description for it's representation
% in the plot you have to be carefull, because sometimes the inner and the
% outher borders are swapped, so if you notice that this is the case you
% have to simply swap the names between, from x_in to x_out and viceversa
% same thing for the y coordinate

w_right = data.w_tr_right_m;
w_left = data.w_tr_left_m;
dx = gradient(x);
dy = gradient(y);
norms = sqrt(dx.^2 + dy.^2);
nx = -dy ./ norms;
ny = dx ./ norms;

% switch here if the borders are inverted:

x_in = (x + nx .* w_left);
y_in = (y + ny .* w_left);


x_out = (x - nx .* w_right);
y_out = (y - ny .* w_right);

%save(path+'track.mat', 'x_in', 'x_out', 'y_in', 'y_out');

% As a last operation, we're extracting the optimal race line from another
% .csv file so that we have always two different paths for our cars on each
% track

data = readtable(path+'raceline.csv');
xp_column = data.x_m;
yp_column = data.y_m;

xp = xp_column';
yp = yp_column';

Pp = [xp_column, yp_column];

%save(path+'optimal_trajectory.mat', 'Pp', 'xp', 'yp');

figure;
hold on;
plot(xp, yp, 'k', 'LineWidth', 2);               % race line
plot(x_out, y_out, 'r--', 'LineWidth', 1.5);     % outher border
plot(x_in, y_in, 'b--', 'LineWidth', 1.5);       % inner border
axis equal;
xlabel('x [m]');
ylabel('y [m]');
title('Circuit');
legend('Race line', 'Outher_border', 'Inner_border');
