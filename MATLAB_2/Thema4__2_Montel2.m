%% Thema 4 ii) me8odo Lyapunov me meikth domh gia online ektimhsh parametrwn

clear;
clc;
close all;


%% pragmatiko systhma

% eisodos u kai parametroi methodou
u = @(t) 1.5 * sin(2 * pi * t) * exp(-3*t);
a = 0.5;
b = 2;

gamma1 = 318;
gamma2 = 89;

theta_m = 5;

% arxikh syn8hkh
x0 = [0, 0, 0, 0]';
t_span = 0:0.01:10;


%% prosomoiwsh kai ektimhsh parametrwn 

% lysh diaforikou systhmatos gia xroniko diasthma t_span
[t, x] = ode15s(@(t, x) system_equationsV6(t, x, u, a, b, gamma1, gamma2, theta_m), t_span, x0);

y = x(:, 1);
y_hat = x(:, 2);
theta_hat1 = x(:, 3);
theta_hat2 = x(:, 4);

% telikes ektimwmenes times parametrwn a, b
a_hat = theta_hat1(:);
b_hat = theta_hat2(:);

a_hat(length(t_span))
b_hat(length(t_span))


%% aksiologhsh basei metrikwn sfalmatos

% error
error = zeros(length(t_span), 1);
for i = 1:length(t_span)
    error(i) = (y(i) - y_hat(i));
end


%% grafikes parastaseis me plots

% kata poso sygklinoyn oi ektimwmenes parametroi tou montelou me to systhma
fig1 = printer_params_est(t_span, a_hat, b_hat, a, b);

% h eksodos y tou systhmatos kai tou montelou
fig2 = printer_outputs(t_span, y, y_hat);
                                                                                              
% sfalma
fig3 = printer_error(t_span, error, true);

saveas(fig1, 'prob4_ii_params_est.png')
saveas(fig2, 'prob4_ii_outputs.png')
saveas(fig3, 'prob4_ii_error.png')