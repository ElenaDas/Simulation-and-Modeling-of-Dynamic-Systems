%% Thema 3 : me8odo Lyapunov me meikth domh gia online ektimhsh parametrwn se deuterhs takshs systhma
clear;
clc;
close all;

%% pragmatiko systhma

% eisodos u kai parametroi methodou
u = @(t) 4 * sin(pi * t) + 2 * sin(8 * pi * t);

% System matrices
a11 = -1;
a12 = 1;
a21 = -4;
a22 = 0;

b1 = 2;
b2 = 1;

A = [a11, a12; a21, a22];
B = [b1; b2];

gamma1 = 15;
gamma2 = 25;
theta_m = [-20, -13; -20, -18];

% arxikh syn8hkh
x0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]';
t_span = 0:0.01:80;


%% prosomoiwsh kai ektimhsh parametrwn 

% lysh diaforikou systhmatos gia xroniko diasthma t_span
[t, x] = ode15s(@(t, x) system_equationsV4(t, x, A, B, gamma1, gamma2, theta_m, u), t_span, x0);

y1 = x(:, 1);
y2 = x(:, 2);
y1_hat = x(:, 9);
y2_hat = x(:, 10);

a11_hat = x(:, 3);
a12_hat = x(:, 4);
a21_hat = x(:, 5);
a22_hat = x(:, 6);

b1_hat = x(:, 7);
b2_hat = x(:, 8);

% telikes ektimwmenes times parametrwn A, B
x(length(t_span), 3:8)

%% aksiologhsh basei metrikwn sfalmatos

% error
error = zeros(length(t_span), 1);
for i = 1:length(t_span)
    error(i) = (y1(i) - y1_hat(i)) + (y2(i) - y2_hat(i)) ;
end

%% grafikes parastaseis me plots

% kata poso sygklinoyn oi ektimwmenes parametroi tou montelou me to systhma
fig0 = printer_params_est_arrayA(t_span, a11_hat, a12_hat, a21_hat, a22_hat, a11, a12, a21, a22);
fig1 = printer_params_estV2(t_span, b1_hat, b2_hat, b1, b2);


% h eksodos y tou systhmatos kai tou montelou
fig2 = printer_outputs(t_span, y1, y1_hat);
fig4 = printer_outputs(t_span, y2, y2_hat);
                          
% meso tetragwniko sfalma
fig3 = printer_error(t_span, error, true);

saveas(fig0, 'prob3_params_est_arrayA.png')
saveas(fig1, 'prob3_params_est.png')
saveas(fig2, 'prob3_outputsY1.png')
saveas(fig4, 'prob3_outputsY2.png')
saveas(fig3, 'prob3_error.png')