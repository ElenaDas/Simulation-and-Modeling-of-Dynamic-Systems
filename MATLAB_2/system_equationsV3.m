function dx = system_equationsV3(t, x, a, b, gamma1, gamma2, theta_m, u, n)
% x = [y, y_hat, theta_hat1, theta_hat2]
    xin = x(1) + n(t);
    e = xin - x(2);
    dx(1) = -a * x(1) + b * u(t);
    dx(2) = -x(3) * x(2) + x(4) * u(t) + theta_m * e;
    dx(3) = -gamma1 * e * xin;
    dx(4) = gamma2 * e * u(t);
    dx = dx';
end