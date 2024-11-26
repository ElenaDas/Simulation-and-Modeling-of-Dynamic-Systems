function fig = printer_params_est_arrayA(t_span, param11, param12, param21, param22, target11, target12, target21, target22)
% param1, param2, target1, target2
    fig = figure;
    plot(t_span, param11, 'Linewidth', 1);
    hold on;
    plot(t_span, param12, 'Linewidth', 1);
    hold on;
    plot(t_span, param21, 'Linewidth', 1);
    hold on;
    plot(t_span, param22, 'Linewidth', 1);
    hold on;
    
    line([t_span(1), t_span(length(t_span))], [target11, target11], 'Color', 'magenta', 'LineStyle','--', 'Linewidth', 1);
    hold on;
    line([t_span(1), t_span(length(t_span))], [target12, target12], 'Color', 'cyan', 'LineStyle','--', 'Linewidth', 1);
    hold on;
    line([t_span(1), t_span(length(t_span))], [target21, target21], 'Color', 'magenta', 'LineStyle','--', 'Linewidth', 1);
    hold on;
    line([t_span(1), t_span(length(t_span))], [target22, target22], 'Color', 'cyan', 'LineStyle','--', 'Linewidth', 1);
    
    legend({'$\hat{a_{11}}$', '$\hat{a_{12}}$', '$\hat{a_{21}}$', '$\hat{a_{22}}$'}, 'Interpreter', 'latex');
    xlabel('$t(sec)$', 'interpreter', 'latex', 'FontWeight', 'bold');
    ylabel('Parameters of array A estimation');
end