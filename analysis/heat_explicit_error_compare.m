function heat_explicit_error_compare()
    clear all; clc; close all;
    alpha = 1.0;
    T = 1.0;   % final time
    fixed_dt = 1e-3;

    nx_list = [11, 21, 41, 81, 161];
    err_fixed = zeros(size(nx_list));
    err_cfl   = zeros(size(nx_list));

    for k = 1:length(nx_list)
        nx = nx_list(k);
        dx = 1 / (nx - 1);
        dt_cfl = 0.4 * dx^2 / alpha;

        % Case A: fixed dt
        [x1, u_num1] = heat_solver_explicit(nx, fixed_dt, T, alpha);
        u_ex1 = exact_solution(x1, T);
        err_fixed(k) = max(abs(u_num1 - u_ex1));

        % Case B: CFL-scaled dt
        [x2, u_num2] = heat_solver_explicit(nx, dt_cfl, T, alpha);
        u_ex2 = exact_solution(x2, T);
        err_cfl(k) = max(abs(u_num2 - u_ex2));
    end

    % Plotting
    figure;
    loglog(nx_list, err_fixed, 'ro-', 'LineWidth', 1.5); hold on;
    loglog(nx_list, err_cfl,   'bs-', 'LineWidth', 1.5);
    grid on; xlabel('Number of spatial points (n_x)');
    ylabel('Max error at T=1');
    legend('Fixed \Deltat', 'CFL-scaled \Deltat', 'Location', 'southwest');
    title('Error vs. Spatial Refinement for Explicit Euler');
end

function [x, u] = heat_solver_explicit(nx, dt, T, alpha)
    dx = 1 / (nx - 1);
    x = linspace(0, 1, nx)';
    nt = round(T / dt);
    r = alpha * dt / dx^2;

    % Stability warning
    if r > 0.5
        warning('r = %.3f exceeds stability limit (0.5).', r);
    end

    u = sin(pi * x);      % Initial condition
    u(1) = 0; u(end) = 0; % Dirichlet BCs

    for n = 1:nt
        u_new = u;
        u_new(2:end-1) = u(2:end-1) + r * (u(3:end) - 2*u(2:end-1) + u(1:end-2));
        u = u_new;
    end
end

function u_exact = exact_solution(x, t)
    u_exact = sin(pi * x) * exp(-pi^2 * t);
end