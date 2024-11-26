%% Κύκλωμα RLC

clear;
clc;
close all;

tic
% Τάσεις εισόδου u1, u2
u1 = @(t) 3 * sin(pi*t);
u2 = 2.5;

% Μηδενικές αρχικές συνθήκες
[x0(1), x0(2)] = v(0);

% Δειγματοληψία για ένα λεπτό ώστε να φανεί σίγουρα του μεταβατικό φαινόμενο 
step = 1e-5;
t_linspace = 0:step:60;
N = length(t_linspace);
y_out = zeros(N, 2);

%% Μετρήσεις αρχείου v.p 
for i = 1:N
    [V_R, V_C] = v(t_linspace(i));
    y_out(i, 1) = V_C;
    y_out(i, 2) = V_R;
end

V_C = y_out(:, 1);

% Τεχνητά σφάλματα για το ερώτημα (β)
% V_C(100000) =  V_C(100000) + 100 *  V_C(100000);
% V_C(200000) =  V_C(200000) + 200 *  V_C(200000);
% V_C(300000) =  V_C(300000) + 300 *  V_C(300000);
% V_C(400000) =  V_C(400000) + 200 *  V_C(400000);
% V_C(500000) =  V_C(500000) + 100 *  V_C(500000);

U1 = double(u1(t_linspace))';
U2 = ones(N, 1) .* u2;


%% Μέθοδος Ελαχίστων Τετραγώνων  
% Πόλοι ευσταθούς φίλτρου                                                   
lambda = [120, 180];

% Διάνυσμα παραμέτρων ζ για το μοντελο δυναμικου συστηματος μετα τη γραμμικη παραμετροποιηση 
z1 = lsim(tf([-1, 0], [1, lambda(1) + lambda(2), lambda(1) * lambda(2)]), V_C, t_linspace');
z2 = lsim(tf(-1, [1, lambda(1) + lambda(2), lambda(1) * lambda(2)]), V_C, t_linspace');
z3 = lsim(tf([1, 0], [1, lambda(1) + lambda(2), lambda(1) * lambda(2)]), U1, t_linspace');
z4 = lsim(tf(1, [1, lambda(1) + lambda(2), lambda(1) * lambda(2)]), U1, t_linspace');
z5 = lsim(tf([1, 0], [1, lambda(1) + lambda(2), lambda(1) * lambda(2)]), U2, t_linspace');
z6 = lsim(tf(1, [1, lambda(1) + lambda(2), lambda(1) * lambda(2)]), U2, t_linspace');
zeta = [z1, z2, z3, z4, z5, z6]';

% Εκτίμηση παραμέτρων συστήματος
theta = optimal_params(y_out, zeta) + [lambda(1) + lambda(2); lambda(1) * lambda(2); 0; 0; 0; 0];
RC_inv = theta(1);
LC_inv = theta(2);


%% ODE 
[t, V_C_estimation] = ode45(@(t, Vc_approx) system_equations(t, Vc_approx, RC_inv, LC_inv, u1, u2), t_linspace, x0);
sum = 0;
for i = 1: N
    sum = sum + abs((V_C(i) - V_C_estimation(i, 1)));
end
V_C_total_error = sum / N;

%V_R = y_out(:, 2); den douleuei gia outliers
V_R = U1 + U2 - V_C;
V_R_estimation = U1 + U2 - V_C_estimation(i, 1);
for i = 1: N
    sum = sum + abs((V_R(i) - V_R_estimation(i)));
end
V_R_total_error = sum / N;
%% Plot V_C
figure(1);
plot(t_linspace, V_C, 'Linewidth', 0.1);
ylabel('Output $y(t) = V_C(t)$', 'interpreter', 'latex');
xlabel('$t(sec)$', 'interpreter', 'latex');

%% Plot εκτιμώμενης V_C
figure(2);
plot(t_linspace, V_C_estimation(:, 1), 'Linewidth', 0.1);
ylabel('Estimated output $\hat{V_C}(t)$', 'interpreter', 'latex');
xlabel('$t(sec)$', 'interpreter', 'latex');

%% Σφάλμα μεταξύ V_C και εκτιμώμενης V_C
VC_error = V_C - V_C_estimation(:, 1);

%% Plot σφάλματος μεταξύ V_C και εκτίμησης V_C
figure(3);
plot(t_linspace, VC_error, 'Linewidth', 0.1);
ylabel('Error $V_C - \hat{V}_C$', 'interpreter', 'latex');
xlabel('$t(sec)$', 'interpreter', 'latex');

%% Υπολογισμός συνολικού σφάλματος V_C και εκτίμησης V_C
VC_total_error = mean(abs(VC_error));
disp(['Total error between actual V_C and estimated V_C: ', num2str(VC_total_error)]);

%% Plot V_R
figure(4);
plot(t_linspace, V_R, 'Linewidth', 0.1);
ylabel('Voltage of resistance $V_R(t)$', 'interpreter', 'latex');
xlabel('$t(sec)$', 'interpreter', 'latex');

%% Plot εκτίμησης V_R
figure(5);
plot(t_linspace, V_R_estimation, 'Linewidth', 0.1);
ylabel('Estimated voltage of resistance $\hat{V_R}(t)$', 'interpreter', 'latex');
xlabel('$t(sec)$', 'interpreter', 'latex');

%% Εκτύπωση παραμέτρων συστήματος 
fprintf('1/RC = %f\n',theta(1) +lambda(1));
fprintf('1/LC = %f',theta(2) + lambda(2));

toc
%% Εξισώσεις του συστήματος 
function dx = system_equations(t, x, RC_inv, LC_inv, u1, u2)
    dx(1)  = x(2);
    dx(2) = -RC_inv * x(2) - LC_inv * x(1) + LC_inv * u2 + RC_inv * u1(t);
    dx = dx';
end

%% Βέλιστες παράμετροι για κάθε στήλη του ζ, MSE 
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

