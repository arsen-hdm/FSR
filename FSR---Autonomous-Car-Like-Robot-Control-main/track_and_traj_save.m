path="Circuito-Monza\";

data = readtable(path+'Monza.csv');
x = data.x_m;
y = data.y_m;


w_right = data.w_tr_right_m;
w_left = data.w_tr_left_m;
dx = gradient(x);
dy = gradient(y);
norms = sqrt(dx.^2 + dy.^2);
nx = -dy ./ norms;
ny = dx ./ norms;

x_out = x + nx .* w_left;
y_out = y + ny .* w_left;


x_in = x - nx .* w_right;
y_in = y - ny .* w_right;

%save('track_matlab.mat', 'x_in', 'x_out', 'y_in', 'y_out');
path="Circuito-Monza\";

data = readtable(path+'Monza_raceline.csv');
xp_colonna = data.x_m;
yp_colonna = data.y_m;

xp = xp_colonna';
yp = yp_colonna';

Pp = [xp_colonna, yp_colonna];

%save('Point_new_map.mat', 'Pp', 'xp', 'yp');

figure;
hold on;
plot(xp, yp, 'k', 'LineWidth', 2);               % Linea centrale
plot(x_out, y_out, 'r--', 'LineWidth', 1.5); % Bordi sinistri
plot(x_in, y_in, 'b--', 'LineWidth', 1.5); % Bordi destri
axis equal;
title('Tracciato del circuito');
legend('Linea centrale', 'Bordo sinistro', 'Bordo destro');
