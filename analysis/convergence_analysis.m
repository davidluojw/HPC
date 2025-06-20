clear all; clc; close all;

h = [1/4, 1/5, 1/10,1/20, 1/30, 1/40, 1/50,1/60, 1/70, 1/80];
log_h = log(h);


k = [1/16000, 1/18000, 1/20000, 1/22000, 1/24000,...
     1/26000, 1/28000, 1/30000, 1/32000, 1/34000];
log_k = log(k);

% =========================================================================
error = [1.540009777446680e-07, 9.289261613851114e-08, 2.224928437705465e-08, 6.368294799133563e-09,3.501335844328912e-09, 2.503289494009767e-09,... 
         2.042279879377198e-09, 1.792104802785557e-09, 1.641342091569537e-09, 1.543525492048502e-09];
log_error = log(error);

error_i = [];
log_error_i = log(error_i);

error_t = [8.908352128169024e-09, 8.065563346262265e-09, 7.388858450700767e-09, 6.834856096857672e-09, 6.372957288599645e-09, ...
           5.981956355395020e-09, 5.646693400789535e-09, 5.356043082744733e-09, 5.101656182359731e-09, 4.877144539060331e-09];
log_error_t = log(error_t);

% =========================================================================


figure;
h1 = plot(log_h, log_error,'b-', 'LineWidth', 2);
xlabel("log-\Deltax");
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

% =========================================================================
figure;
h2 = plot(log_k, log_error_t,'b-', 'LineWidth', 2);
xlabel("log-\Deltat");
ylabel("log-error");
legend('error','Location', 'Best', 'FontSize', 14, 'Box', 'on');

rate_error_t = zeros(length(log_k)-1,1);

for ii = 1:length(log_k)-1
    rate_error_t(ii) = (log_error_t(ii+1) - log_error_t(ii))/(log_k(ii+1)-log_k(ii));
end

% add slope text
for seg = 1:length(log_h)-1
    x_mid = mean(log_k(seg:seg+1));
    y_mid_t = mean(log_error_t(seg:seg+1)) + 0.05;
    text(x_mid, y_mid_t, sprintf('Slope: %.4f',rate_error_t(seg)),...
        'Color', 'b',...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'bottom',...
        'FontSize', 12);
end
legend('show'); 
hold off;

% =========================================================================


figure;
h1 = plot(log_h, log_error_i,'b-', 'LineWidth', 2);
xlabel("log-\Deltax");
ylabel("log-error");
legend('error','Location', 'Best', 'FontSize', 14, 'Box', 'on');

% rate of error
rate_error_i = zeros(length(log_h)-1,1);
for ii = 1:length(log_h)-1
    rate_error_i(ii) = (log_error_i(ii+1) - log_error_i(ii))/(log_h(ii+1)-log_h(ii));
end


% add slope text
for seg = 1:length(log_h)-1
    x_mid = mean(log_h(seg:seg+1));
    y_mid = mean(log_error_i(seg:seg+1)) + 0.2;
    text(x_mid, y_mid, sprintf('Slope: %.4f',rate_error_i(seg)),...
        'Color', 'b',...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'bottom',...
        'FontSize', 12);
end
legend('show'); 
hold off;

% =========================================================================

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
