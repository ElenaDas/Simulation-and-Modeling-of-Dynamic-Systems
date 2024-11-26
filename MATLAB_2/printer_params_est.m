function fig = printer_params_est(t_span, param1, param2, target1, target2)
    fig = figure;
    plot(t_span, param1, 'Linewidth', 1);
    hold on;
    plot(t_span, param2, 'Linewidth', 1);
    hold on;
    
    line([t_span(1), t_span(length(t_span))], [target1, target1], 'Color', 'magenta', 'LineStyle','--', 'Linewidth', 1);
    hold on;
    line([t_span(1), t_span(length(t_span))], [target2, target2], 'Color', 'cyan', 'LineStyle','--', 'Linewidth', 1);
    
    legend({'$\hat{a}$', '$\hat{b}$', '$a$', '$b$'}, 'Interpreter', 'latex');
    xlabel('$t(sec)$', 'interpreter', 'latex', 'FontWeight', 'bold');
    ylabel('$a, b, \hat{a}(t), \hat{b}(t)$', 'interpreter', 'latex', ...
                                 'FontWeight','bold', 'FontSize', 12);
end