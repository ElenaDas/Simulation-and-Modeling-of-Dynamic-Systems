function fig = printer_outputs(t_span, y, y_hat)
    fig = figure;
    plot(t_span, y, 'Linewidth', 1);
    hold on;
    plot(t_span, y_hat, 'Linewidth', 1);
    legend({'$x$', '$\hat{x}$'}, 'Interpreter', 'latex');
    xlabel('$t(sec)$', 'interpreter', 'latex', 'FontWeight', 'bold');
    ylabel('$x(t), \hat{x}(t)$', 'interpreter', 'latex', 'FontWeight', 'bold');
    
end