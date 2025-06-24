clear all; clc; close all;

h = [1/20, 1/40, 1/60, 1/80, 1/100];
log_h = log(h);


k = [1/50, 1/500, 1/5000, 1/50000];
log_k = log(k);


% =========================================================================
h = [1/20, 1/40, 1/60, 1/80, 1/100];
log_h = log(h);

error_ex_spatial = [ 0.0027879554355500,  0.0006966557442398, 0.0003096218467706, 0.0001741779746242, 0.0001114899679588];
log_error_ex_sp = log(error_ex_spatial);

figure;
h1 = plot(log_h, log_error_ex_sp,'b-', 'LineWidth', 2);
xlabel("log-\Deltax");
ylabel("log-error");
legend('error','Location', 'Best', 'FontSize', 14, 'Box', 'on');

% rate of error
rate_error_ex_sp = zeros(length(log_h)-1,1);
for ii = 1:length(log_h)-1
    rate_error_ex_sp(ii) = (log_error_ex_sp(ii+1) - log_error_ex_sp(ii))/(log_h(ii+1)-log_h(ii));
end


% add slope text
for seg = 1:length(log_h)-1
    x_mid = mean(log_h(seg:seg+1));
    y_mid = mean(log_error_ex_sp(seg:seg+1)) + 0.2;
    text(x_mid, y_mid, sprintf('%.4f',rate_error_ex_sp(seg)),...
        'Color', 'b',...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'bottom',...
        'FontSize', 12);
end
legend('show'); 
hold off;

% =========================================================================
h = [1/20, 1/40, 1/60, 1/80, 1/100];
log_h = log(h);

error_im_spatial = [ 0.0027871483886043, 0.0006958398002474, 0.0003088042278093, 0.0001733597864412, 0.0001106715025248];
log_error_im_sp = log(error_im_spatial);  


figure;
h1 = plot(log_h, log_error_im_sp,'b-', 'LineWidth', 2);
xlabel("log-\Deltax");
ylabel("log-error");
legend('error','Location', 'Best', 'FontSize', 14, 'Box', 'on');

% rate of error
rate_error_im_sp = zeros(length(log_h)-1,1);
for ii = 1:length(log_h)-1
    rate_error_im_sp(ii) = (log_error_im_sp(ii+1) - log_error_im_sp(ii))/(log_h(ii+1)-log_h(ii));
end


% add slope text
for seg = 1:length(log_h)-1
    x_mid = mean(log_h(seg:seg+1));
    y_mid = mean(log_error_im_sp(seg:seg+1)) + 0.2;
    text(x_mid, y_mid, sprintf('%.4f',rate_error_im_sp(seg)),...
        'Color', 'b',...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'bottom',...
        'FontSize', 12);
end
legend('show'); 
hold off;


% =========================================================================
k = [1/10000, 1/20000, 1/40000, 1/80000];
log_k = log(k);

error_ex_temporal = [0.0000020121666123, 0.0000010060695050, 0.0000005030313020, 0.0000002515147885];
log_error_ex_tp = log(error_ex_temporal);

figure;
h2 = plot(log_k, log_error_ex_tp,'b-', 'LineWidth', 2);
xlabel("log-\Deltat");
ylabel("log-error");
legend('error','Location', 'Best', 'FontSize', 14, 'Box', 'on');

rate_error_ex_tp = zeros(length(log_k)-1,1);

for ii = 1:length(log_k)-1
    rate_error_ex_tp(ii) = (log_error_ex_tp(ii+1) - log_error_ex_tp(ii))/(log_k(ii+1)-log_k(ii));
end

% add slope text
for seg = 1:length(log_k)-1
    x_mid = mean(log_k(seg:seg+1));
    y_mid_t = mean(log_error_ex_tp(seg:seg+1)) + 0.05;
    text(x_mid, y_mid_t, sprintf('%.4f',rate_error_ex_tp(seg)),...
        'Color', 'b',...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'bottom',...
        'FontSize', 12);
end
legend('show'); 
hold off;

% =========================================================================
k = [1/10000, 1/20000, 1/40000, 1/80000];
log_k = log(k);

error_im_temporal = [0.0000020120559543, 0.0000010060412947, 0.0000005030242415, 0.0000002515127275];
log_error_im_tp = log(error_im_temporal);


figure;
h2 = plot(log_k, log_error_im_tp,'b-', 'LineWidth', 2);
xlabel("log-\Deltat");
ylabel("log-error");
legend('error','Location', 'Best', 'FontSize', 14, 'Box', 'on');

% rate of error
rate_error_im_tp = zeros(length(log_k)-1,1);
for ii = 1:length(log_k)-1
    rate_error_im_tp(ii) = (log_error_im_tp(ii+1) - log_error_im_tp(ii))/(log_k(ii+1)-log_k(ii));
end


% add slope text
for seg = 1:length(log_k)-1
    x_mid = mean(log_k(seg:seg+1));
    y_mid = mean(log_error_im_tp(seg:seg+1)) + 0.2;
    text(x_mid, y_mid, sprintf('%.4f',rate_error_im_tp(seg)),...
        'Color', 'b',...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'bottom',...
        'FontSize', 12);
end
legend('show'); 
hold off;

% =========================================================================