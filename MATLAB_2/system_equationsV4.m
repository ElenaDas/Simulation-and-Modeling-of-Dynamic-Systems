function dx = system_equationsV4(t, x, A, B, gamma1, gamma2, theta_m, u)
% x = [y1, y2, a11_hat, a12_hat, a21_hat, a22_hat, b1_hat, b2_hat, y1_hat, y2_hat]
    e(1) = x(1) - x(9);
    e(2) = x(2) - x(10);
    dx(1) = A(1, 1) * x(1) + A(1, 2) * x(2) + B(1) * u(t);
    dx(2) = A(2, 1) * x(1) + A(2, 2) * x(2) + B(2) * u(t);
    dx(3) = gamma1 * e(1) * x(1);
    dx(4) = gamma1 * e(1) * x(2);
    dx(5) = gamma1 * e(2) * x(1);
    dx(6) = gamma1 * e(2) * x(2);
    dx(7) = gamma2 * e(1) * u(t);
    dx(8) = gamma2 * e(2) * u(t);
    dx(9) = x(3) * x(9) + x(4) * x(10) + x(7) * u(t) - (theta_m(1, 1) * e(1) + theta_m(1, 2) * e(2));
    dx(10) = x(5) * x(9) + x(6) * x(10) + x(8) * u(t) - (theta_m(2, 1) * e(1) + theta_m(2, 2) * e(2));
    dx = dx';
end