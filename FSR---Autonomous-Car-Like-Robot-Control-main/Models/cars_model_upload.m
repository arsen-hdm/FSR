l = 2.2;          %m distance between the front and rear wheels
phi_lim = 0.45;   %rad maximum steering angle
vlim = 36;        %m/s maximum heading velocity
wheelR = 0.26;    %m radius of the wheels
wlim = 1.0;       %rad/s maximum steering velocity

save('Fiat_600_model.mat', 'l', 'phi_lim','vlim','wheelR','wlim');