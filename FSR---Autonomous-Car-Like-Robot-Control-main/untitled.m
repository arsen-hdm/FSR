%Almost non linear bicycle:
clc; clear; close all;

% === Parametri simulazione ===
Ts = 0.025;
T_sim = 60;
L = 2.5;         % distanza tra assi
Kp = [1.5, 3.5, 2.0];  % controlli: [Kx, Ky, Ktheta]
tau_delta = 0.3;       % tempo attuazione sterzo (per derivata δ)

% === Caricamento dati ===
load('circuit_Suzuka/optimal_trajectory.mat');   % xp, yp
load('circuit_Suzuka/track.mat');    % x_in, y_in, x_out, y_out

% === Raceline regolare ===
d = sqrt(diff(xp).^2 + diff(yp).^2);
s_track = [0, cumsum(d)];
s_total = s_track(end);
s_uniform = linspace(0, s_total, 5 * floor(T_sim/Ts));

xp_ref = interp1(s_track, xp, s_uniform, 'spline');
yp_ref = interp1(s_track, yp, s_uniform, 'spline');

dx = gradient(xp_ref);
dy = gradient(yp_ref);
ddx = gradient(dx);
ddy = gradient(dy);
theta_ref = unwrap(atan2(dy, dx));

curvature = abs(dx .* ddy - dy .* ddx) ./ (dx.^2 + dy.^2).^(3/2);
curvature(isnan(curvature)) = 0;

% === Stato iniziale (modello bicycle: x, y, θ, δ) ===
x = xp_ref(1) + 0.5;
y = yp_ref(1) + 0.5;
theta = theta_ref(1) + deg2rad(20);
delta = 0;  % angolo di sterzo

traj = zeros(4, length(s_uniform));
traj(:,1) = [x; y; theta; delta];
v_old = 0; delta_dot_old = 0;

% === Simulazione ===
for k = 1:length(s_uniform)-1
    x_d = xp_ref(k); y_d = yp_ref(k); theta_d = theta_ref(k);
    dx = xp_ref(k+1) - xp_ref(k);
    dy = yp_ref(k+1) - yp_ref(k);

    curv_k = curvature(k);
    v_d = max(2.5, min(8, 6.5 / (1 + 20 * curv_k)));
    omega_d = (theta_ref(k+1) - theta_ref(k)) / Ts;

    dx_e = x_d - x;
    dy_e = y_d - y;
    ex = cos(theta)*dx_e + sin(theta)*dy_e;
    ey = -sin(theta)*dx_e + cos(theta)*dy_e;
    etheta = wrapToPi(theta_d - theta);

    % === Controllo Almost-Nonlinear adattato a bicycle ===
    v = v_d * cos(etheta) + Kp(1)*ex;
    omega = omega_d + Kp(2)*ey + Kp(3)*sin(etheta);

    % Calcolo δ_des da ω
    delta_des = atan(L * omega / max(v, 0.1));

    % Modello attuatore sterzo (1° ordine)
    delta_dot = (delta_des - delta) / tau_delta;
    delta = delta + delta_dot * Ts;

    % Dinamica bicycle
    x = x + v * cos(theta) * Ts;
    y = y + v * sin(theta) * Ts;
    theta = theta + (v / L) * tan(delta) * Ts;

    traj(:,k+1) = [x; y; theta; delta];
    v_old = v; delta_dot_old = delta_dot;
end

% === Plot finale ===
figure;
hold on; axis equal; grid on;
plot(xp_ref, yp_ref, 'k--', 'LineWidth', 2, 'DisplayName', 'Raceline');
plot(traj(1,:), traj(2,:), 'b', 'LineWidth', 2, 'DisplayName', 'Bicycle Trajectory');
plot(x_in, y_in, 'r--', 'LineWidth', 1, 'DisplayName','Bordo interno');
plot(x_out, y_out, 'r--', 'LineWidth', 1, 'DisplayName','Bordo esterno');
legend('Location', 'best');
title('Bicycle Model Tracking della Raceline');
xlabel('X [m]'); ylabel('Y [m]');