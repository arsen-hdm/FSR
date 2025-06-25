% Punti della traiettoria (esempio)
x = xp;  % vettore delle coordinate x
y = yp;  % vettore delle coordinate y

% Parametri veicolo
v_max = 36;                  % m/s
delta_max = deg2rad(30);     % angolo di sterzata massimo [rad]
omega_delta_max = 0.5;       % velocità massima di sterzata [rad/s]
a_lat_max = 3;               % accelerazione laterale massima [m/s^2]
a_long_max = 2.0;            % accelerazione longitudinale massima [m/s^2]
L = 2.2;                     % passo del veicolo [m]

% Calcolo curvilinear abscissa s
dx = diff(x);
dy = diff(y);
ds = sqrt(dx.^2 + dy.^2);
s = [0 cumsum(ds)];

% Derivate numeriche
% Derivate numeriche corrette e compatibili
theta = atan2(dy, dx);           % lunghezza N-1
dtheta = diff(unwrap(theta));    % N-2, uso unwrap per continuità angoli
curvature = [0, dtheta ./ ds(2:end)];  % N-1: 0 iniziale, derivata altrove
kappa = [curvature, 0];          % match dimensione con s


% Rimozione artefatti di curvatura
kappa = medfilt1(kappa, 5);

% Velocità massima per accelerazione laterale
v_lat_max = sqrt(a_lat_max ./ (abs(kappa) + 1e-6));

% Velocità massima per angolo di sterzata
delta_curve = atan(L * kappa);
v_delta_max = v_max * ones(size(delta_curve));
v_delta_max(abs(delta_curve) > delta_max) = ...
    sqrt(abs(tan(delta_max) ./ (abs(kappa(abs(delta_curve) > delta_max)) * L)));

% Velocità massima per velocità di sterzata (grezza)
omega_delta = [0 diff(delta_curve)./ds];
v_omega_delta_max = omega_delta_max ./ (abs(omega_delta) + 1e-6);
v_omega_delta_max = min(v_omega_delta_max, v_max);

% Velocità desiderata finale (limitata da tutte le condizioni)
v_des_raw = min([v_max * ones(size(v_lat_max)); 
                 v_lat_max; 
                 v_delta_max; 
                 v_omega_delta_max], [], 1);

% Smooth della velocità per continuità
v_des_smooth = smoothdata(v_des_raw, 'gaussian', 10);

% Costruzione velocità crescente realistica da zero
v_profile = zeros(size(v_des_smooth));
v_profile(1) = 0;

for i = 2:length(v_profile)
    ds_i = s(i) - s(i-1);
    v_prev = v_profile(i-1);
    v_accel_limit = sqrt(v_prev^2 + 2 * a_long_max * ds_i);
    v_profile(i) = min(v_accel_limit, v_des_smooth(i));
end

% % (Facoltativo) Rallentamento alla fine
% for i = length(v_profile)-1:-1:1
%     ds_i = s(i+1) - s(i);
%     v_next = v_profile(i+1);
%     v_decel_limit = sqrt(v_next^2 + 2 * a_long_max * ds_i);
%     v_profile(i) = min(v_profile(i), v_decel_limit);
% end

% Tempo cumulativo
t = zeros(size(s));
for i = 2:length(s)
    v_avg = (v_profile(i) + v_profile(i-1)) / 2;
    t(i) = t(i-1) + ds(i-1) / max(v_avg, 1e-3);
end

% Plot
figure;
subplot(3,1,1)
plot(s, v_profile); grid on;
xlabel('Curvilinear abscissa [m]'); ylabel('v_{des} [m/s]');
title('Velocità desiderata realistica');

subplot(3,1,2)
plot(s, t); grid on;
xlabel('Curvilinear abscissa [m]'); ylabel('Tempo [s]');
title('Tempo cumulato');

subplot(3,1,3)
plot(x, y, 'k'); axis equal; grid on;
xlabel('x [m]'); ylabel('y [m]');
title('Traiettoria');

% Sterzata desiderata
delta_des = atan(L * kappa);

% Plot aggiuntivo: Sterzata desiderata
figure;
plot(s, rad2deg(delta_des)); grid on;
xlabel('Curvilinear abscissa [m]');
ylabel('\delta_{des} [deg]');
title('Sterzata desiderata lungo la traiettoria');

theta = atan2(dy, dx);           % orientamento lungo traiettoria (N-1)

% Calcolo velocità angolare orientamento (theta_dot)
dtheta = diff(unwrap(theta)) ./ ds;    % N-2 valori di derivata rispetto a s

% Allineo la dimensione aggiungendo un valore finale (ad esempio ripetendo l'ultimo)
theta_dot = [dtheta, dtheta(end)];     % velocità angolare orientamento

% Calcolo velocità angolare di sterzata (omega = delta_dot)
            % angolo sterzata desiderato
ddelta = diff(delta_des) ./ ds;         % derivata angolo sterzata rispetto a s
omega_ref = [ddelta, ddelta(end)];      % velocità angolare sterzata

% Nota: la derivata è rispetto a s, per passare a derivata rispetto a tempo t
% si usa chain rule: d/dt = d/ds * ds/dt = derivata_s * v_profile

theta_dot_t = theta_dot .* v_profile;   % [rad/s]
omega_ref_t = omega_ref .* v_profile;   % [rad/s]

% Orientamento desiderato (theta)
theta_ref = theta;

% Velocità angolare orientamento desiderata (theta dot)
theta_dot_ref = theta_dot_t;

% Velocità angolare sterzata desiderata (omega)
omega_ref = omega_ref_t;

% Plot
figure;
subplot(5,1,1)
plot(t, v_profile, 'LineWidth', 1.5);
grid on; xlabel('Tempo [s]'); ylabel('Velocità [m/s]');
title('Velocità desiderata (v)');

subplot(5,1,2)
plot(t, delta_ref, 'LineWidth', 1.5);
grid on; xlabel('Tempo [s]'); ylabel('Angolo sterzata \delta [rad]');
title('Angolo di sterzata desiderato (\delta)');

subplot(5,1,3)
plot(t, theta, 'LineWidth', 1.5);
grid on; xlabel('Tempo [s]'); ylabel('\theta [rad]');
title('Orientamento desiderato (\theta)');

subplot(5,1,4)
plot(t, theta_dot, 'LineWidth', 1.5);
grid on; xlabel('Tempo [s]'); ylabel('Velocità angolare \theta\_dot [rad/s]');
title('Velocità angolare orientamento desiderata (\theta dot)');

subplot(5,1,5)
plot(t, omega_ref, 'LineWidth', 1.5);
grid on; xlabel('Tempo [s]'); ylabel('Velocità angolare sterzata \omega [rad/s]');
title('Velocità angolare sterzata desiderata (\omega)');


