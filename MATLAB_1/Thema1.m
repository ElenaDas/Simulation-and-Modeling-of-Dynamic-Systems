%% Κίνηση μάζας με εξωτερικές δυνάμεις, ελατήριο και αποσβεστήρα

clear;
clc;
close all;

tic

m = 8.5; % [kg]
b = 0.65; % [kg/s]
k = 2; % [kg/s^2]
u = @(t) 10 * cos(0.5*pi*t) + 3; % [N] 

% Μηδενικές αρχικές συνθήκες
x0(1) = 0;
x0(2) = 0;

% Δειγματοληψία ανά 0.1 δευτερόλεπτο για διάρκεια 10 δευτερόλεπτα
step = 0.1;
t_linspace = 0:step:10;


%% ODE 
[t, y_out] = ode15s(@(t, y_out)system_equations(t, y_out, m, b, k, u), t_linspace, x0);                                       
U = u(t(:));
Y = y_out(:, 1);

%% Μέθοδος ελαχίστων τετραγώνων(LSM) 
% Πόλοι ευσταθπυς φίλτρου 
lambda = [0.5, 0.5];

% Διάνυσμα παραμέτρων μετά την γραμμική παραμετροποίηση 
syms m_ b_ k_; 
theta_lambda = [b_ / m_ - (lambda(1) + lambda(2)); k_ / m_ - (lambda(1) * lambda(2)); 1 / m_];

% Διάνυσμα εισόδου
% lsim(sys,u,t)
z1 = lsim(tf([-1, 0], [1, lambda(1) + lambda(2), lambda(1) * lambda(2)]), y_out(:, 1), t);
z2 = lsim(tf(-1, [1, lambda(1) + lambda(2), lambda(1) * lambda(2)]), y_out(:, 1), t);
z3 = lsim(tf(1, [1, lambda(1) + lambda(2), lambda(1) * lambda(2)]), U, t);                                        
zeta = [z1, z2, z3]';

% Παράμετροι
op = optimal_params(y_out, zeta);
eqns = (theta_lambda == op);
solved_params = solve(eqns, [m_, b_, k_]);
m_ = double(solved_params.m_);
b_ = double(solved_params.b_);
k_ = double(solved_params.k_);


estimated_params = [m_, b_, k_];
disp('Estimated parameters:');
disp(estimated_params);


fprintf("Least squares method calculated the following parameters of the model: m = %d, b = %d, k = %d\n", m_, b_, k_);
%% ODE 
[t, y_out_estimation] = ode15s(@(t, y_out_estimation) system_equations(t, y_out_estimation, m_, b_, k_, u), t_linspace, x0);
error = zeros(length(y_out), 1); 
N = length(y_out);
for i = 1:N
    error(i) = abs((y_out(i) - y_out_estimation(i)));
end

%% Plots εισόδου, εξόδου και σφάλματος μεταξύ τους
figure(1);
plot(t, y_out(:, 1), 'Linewidth', 1);
ylabel('$y(t)$', 'Interpreter', 'latex');
xlabel('$t(sec)$', 'Interpreter', 'latex');
title('Output')

figure(2);
plot(t, y_out_estimation(:, 1), 'Linewidth', 1);
ylabel('$\hat{y}(t)$', 'Interpreter', 'latex');
xlabel('$t(sec)$', 'Interpreter', 'latex');
title('Estimated Output');

figure(3);
plot(t, error, 'Linewidth', 1);
ylabel('e = $y(t) - \hat{y}(t)$','Interpreter', 'latex');
xlabel('$t(sec)$', 'Interpreter', 'latex');
title('Deviation between actual and estimated Output')

toc
%% Εξισώσεις του συστήματος
function dx = system_equations(t, x, m, b, k, u)
    dx(1) = x(2);
    dx(2) = (1 / m) * (-k * x(1) - b * x(2) + u(t));
    dx = dx';
end

%% Βέλτιστες παράμετροι για κάθε στήλη του ζ, MSE 
function theta = optimal_params(y_out, zeta)
    sum_denominator = 0;
    sum_numerator = 0;
    N = length(y_out);
    for i = 1:N
        sum_denominator = sum_denominator + zeta(:, i) * zeta(:, i)';
        sum_numerator = sum_numerator + zeta(:, i) * y_out(i, 1);
    end
    
    theta = sum_denominator \ sum_numerator;
end