function fig = printer_params_estV2(t_span, param1, param2, target1, target2)
    fig = figure;
    plot(t_span, param1, 'Linewidth', 1);
    hold on;
    plot(t_span, param2, 'Linewidth', 1);
    hold on;
    
    line([t_span(1), t_span(length(t_span))], [target1, target1], 'Color', 'magenta', 'LineStyle','--', 'Linewidth', 1);
    hold on;
    line([t_span(1), t_span(length(t_span))], [target2, target2], 'Color', 'cyan', 'LineStyle','--', 'Linewidth', 1);
    
    legend({'$\hat{b_1}$', '$\hat{b_2}$', '$b_1$', '$b_2$'}, 'Interpreter', 'latex');
    xlabel('$t(sec)$', 'interpreter', 'latex', 'FontWeight', 'bold');
    ylabel('Parameter estimation');
end