function dx = system_equationsV5(t, x, u, a, b, gamma1, gamma2, theta_m)
% x = [y, y_hat, theta_hat1, theta_hat2]
   
    e = x(1) - x(2);
    dx(1) = -a * 0.5 * sin(x(1)) * x(1)   + b * u(t);
    dx(2) = -x(3) *  0.5 * sin(x(1)) * x(1)  + x(4) * u(t) + theta_m * (x(1)-x(2));
    dx(3) = - gamma1 * e * 0.5 * sin(x(1)) * x(1);
    dx(4) = gamma2 * e * u(t);
    dx = dx';
end