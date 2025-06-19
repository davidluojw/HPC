clear all; clc; close all;

figure;
h = [1/5, 1/10,1/20, 1/30, 1/40, 1/50,1/60, 1/70, 1/80];
log_h = log(h);



% =========================================================================
error = [9.289261613851114e-08, 2.224928437705465e-08, 6.368294799133563e-09,3.501335844328912e-09, 2.503289494009767e-09,... 
         2.042279879377198e-09, 1.792104802785557e-09, 1.641342091569537e-09, 1.543525492048502e-09];
log_error = log(error);

% =========================================================================



h1 = plot(log_h, log_error,'b-', 'LineWidth', 2);
xlabel("log-h");
ylabel("log-error");
legend('error','Location', 'Best', 'FontSize', 14, 'Box', 'on');

% rate of error
rate_error = zeros(length(log_h)-1,1);
for ii = 1:length(log_h)-1
    rate_error(ii) = (log_error(ii+1) - log_error(ii))/(log_h(ii+1)-log_h(ii));
end

% add slope text
for seg = 1:length(log_h)-1
    x_mid = mean(log_h(seg:seg+1));
    y_mid = mean(log_error(seg:seg+1)) + 0.2;
    text(x_mid, y_mid, sprintf('Slope: %.4f',rate_error(seg)),...
        'Color', 'b',...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'bottom',...
        'FontSize', 12);
end
legend('show'); 
hold off;

fprintf('\n                           Convergence\n');
fprintf('-------------------------------------------------------------------------------\n');
fprintf('\tgrid\te\tlog(e)\n');

for i = 1:length(log_h)
    e = error(i);
    loge = log_error(i); 
    
    % print
    line = sprintf('\t%d\t%.8f\t%.8f\n', i, e, loge);
    fprintf('%s\n', line);
end
