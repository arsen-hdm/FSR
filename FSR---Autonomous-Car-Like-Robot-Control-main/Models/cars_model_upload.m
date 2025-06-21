clc; clear all; close all;

l = 2.2;          %m distance between the front and rear wheels
phi_lim = deg2rad(somethig);   %rad maximum steering angle
vlim = 36;        %m/s maximum heading velocity
wheelR = 0.26;    %m radius of the wheels
wlim = 1;       %rad/s maximum steering velocity

% change the name of the saved .mat file to distinguish between different
% models

save('Fiat_600.mat', 'l', 'phi_lim','vlim','wheelR','wlim');