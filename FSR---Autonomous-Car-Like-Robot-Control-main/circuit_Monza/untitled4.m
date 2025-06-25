% --- INPUT: traiettoria e parametri veicolo ---
x = xp;     % coordinate x della traiettoria
y = yp;     % coordinate y della traiettoria

v_max = 36;
delta_max = deg2rad(30);
omega_delta_max = 0.5;
a_lat_max = 3.0;
a_long_max = 2.0;
L = 2.2;

% --- 1. Curvilinear abscissa s ---
dx = diff(x);
dy = diff(y);
ds = sqrt(dx.^2 + dy.^2);
s = [0, cumsum(ds)];

% --- 2. Derivate numeriche ---
theta = atan2(dy, dx);
dtheta = diff(unwrap(theta));
curvature = [0, dtheta ./ ds(2:end)];
kappa = [curvature, 0];
kappa = medfilt1(kappa, 5);

% --- 3. Velocità ammissibili ---
v_lat_max = sqrt(a_lat_max ./ (abs(kappa) + 1e-6));

delta_curve = atan(L * kappa);
v_delta_max = v_max * ones(size(delta_curve));
mask_delta = abs(delta_curve) > delta_max;
v_delta_max(mask_delta) = sqrt(abs(tan(delta_max) ./ (abs(kappa(mask_delta)) * L)));

omega_delta = [0, diff(delta_curve) ./ ds];
v_omega_delta_max = omega_delta_max ./ (abs(omega_delta) + 1e-6);
v_omega_delta_max = min(v_omega_delta_max, v_max);

% --- 4. Profilo velocità desiderato lungo s ---
v_des_raw = min([v_max * ones(size(v_lat_max)); 
                 v_lat_max; 
                 v_delta_max; 
                 v_omega_delta_max], [], 1);
v_tilde = smoothdata(v_des_raw, 'gaussian', 10);  % velocità geometrica

% --- 5. Simula un profilo temporale v_t(t) (es. crescente realistica) ---
v_profile = zeros(size(v_tilde));
v_profile(1) = 0;
for i = 2:length(v_profile)
    ds_i = s(i) - s(i-1);
    v_prev = v_profile(i-1);
    v_accel_limit = sqrt(v_prev^2 + 2 * a_long_max * ds_i);
    v_profile(i) = min(v_accel_limit, v_tilde(i));
end
v_t = v_profile;

% --- 6. Calcolo s_dot(t) = v_t(t) / v_tilde(s) ---
n = length(s);
t = zeros(1, n);
for i = 2:n
    v_avg = (v_t(i) + v_t(i-1)) / 2;
    ds_i = s(i) - s(i-1);
    t(i) = t(i-1) + ds_i / max(v_avg, 1e-3);
end
t = t - t(1);  % assicura che t parta da 0

% Interpola s(t) da t
s_t = @(tq) interp1(t, s, tq, 'pchip');
s_values = s_t(t);

% --- 7. Ricostruzione delle variabili temporali ---
x_t_reference = @(tq) interp1(s, x, s_t(tq), 'spline');
y_t_reference = @(tq) interp1(s, y, s_t(tq), 'spline');

x_values = arrayfun(x_t_reference, t);
y_values = arrayfun(y_t_reference, t);

% Derivata numerica x(t), y(t)
dt = gradient(t);
x_dot = gradient(x_values) ./ dt;
y_dot = gradient(y_values) ./ dt;
v_t_check = sqrt(x_dot.^2 + y_dot.^2);  % controllo v_t

% --- 8. Output come timeseries ---
x_ref = timeseries(x_values, t);
y_ref = timeseries(y_values, t);
v_ref = timeseries(v_t, t);

% --- 9. (facoltativo) Plotta i risultati ---
figure;
subplot(3,1,1); plot(s, v_tilde, 'k', s, v_t, 'b'); legend('v_{tilde}(s)', 'v_t'); title('Velocità');
subplot(3,1,2); plot(t, s_values); title('s(t)');
subplot(3,1,3); plot(x, y, 'k--'); hold on; plot(x_values, y_values, 'b'); title('Traiettoria');

% --- 10. Calcolo theta_d(t) = atan2(y_dot, x_dot) ---
theta_d = unwrap(atan2(y_dot, x_dot));
theta_dot_d = gradient(theta_d) ./ dt;


% --- 11. Calcolo φ_d(t) = atan(L * θ_dot / v_t) ---
phi_d = atan(L * theta_dot_d ./ (v_t + 1e-6));
phi_dot_d = gradient(phi_d) ./ dt;

% --- 12. Output finale come timeseries ---
x_t_reference       = timeseries(x_values, t);
y_t_reference      = timeseries(y_values, t);
x_dot_t_reference   = timeseries(x_dot, t);
y_dot_t_reference   = timeseries(y_dot, t);
theta_t_reference   = timeseries(theta_d, t);
theta_t_dot_reference = timeseries(theta_dot_d, t);
phi_t_reference     = timeseries(phi_d, t);
w_t_reference   = timeseries(phi_dot_d, t);
v_t_ts       = timeseries(v_t, t);

theta_values = arrayfun(@(t) theta_t_reference(t),t);

% --- 13. Visualizzazione finale (facoltativa) ---
figure;
subplot(3,2,1); plot(x_values); title('x_d(t)');
subplot(3,2,2); plot(y_values); title('y_d(t)');
subplot(3,2,3); plot(t, x_dot); title('x\_dot\_d(t)');
subplot(3,2,4); plot(t, y_dot); title('y\_dot\_d(t)');
subplot(3,2,5); plot(theta_t_reference); title('\theta_d(t)');
subplot(3,2,6); plot(phi_t_reference); title('\phi_d(t)');
tf = 300


