%% 1.1

clc
close
clear
% Setting parameters for all the problems
g = 1; % gravity acceleration
x_initial = 0;
x_final = 2;
T = 0.5; % Final time
CFL = 0.25; % CFL constant
f = @(q) [q(2,:) + 0*q(1,:); q(2,:).^2./q(1,:) + 0.5*g*q(1,:).^2]; % flux
u = 0.25; %velocity
% Setting initial conditions
h0 = @(x) 1 + 0.5*sin(pi*x);
m0 = @(x) u*h0(x);

% Defining source term
S = @(x,t) [pi/2 * (u-1)*cos(pi*(x-t));...
            pi/2 * cos(pi*(x-t)) .* (-u + u^2 + g*h0(x-t))];
% Exact solutions for error analysis
h_ex = @(x,t) h0(x-t);
m_ex = @(x,t) u*h_ex(x,t);
% type of boundary conditions
bc="periodic";
% Defining number of cells for each numerical solution computation
N_vec = 2.^[7,8,9,10];
err_table = zeros(3, length(N_vec));

figure(1)
hold on
% Setting colors for different plots
colors = lines(length(N_vec)+1);

% Solving numerical problem for all different timesteps
for i = 1:length(N_vec)
    N=N_vec(i);
   
    % Solve for Lax-Friedrichs conservative flux method
    [xx, t, q] = LaxFriedrichs(x_initial, x_final, N, T, bc, h0, m0, f, S, CFL, g);
    dx= xx(2)-xx(1);

    % Plotting numerical results for h and m
    subplot(1,2,1)
    plot(xx, q(1,:), "Color", colors(i, :));
    hold on
    subplot(1,2,2)
    plot(xx, q(2,:), "Color", colors(i, :));
    hold on

    % Computing errors for h and m
    m_err_inf = norm(q(2,:)-m_ex(xx, t), "inf");
    h_err_inf = norm(q(1,:)-h_ex(xx, t), "inf");
    err_table(3,i)=m_err_inf;
    err_table(2,i)=h_err_inf;
    err_table(1,i)=dx;
end

% Plotting exact solution for h
subplot(1,2,1)
plot(xx, h_ex(xx, T), 'LineWidth', 1.5, "Color", colors(end, :));
legend_strings = arrayfun(@(x) sprintf('N=%d', x), N_vec, 'UniformOutput', false);
legend([legend_strings, 'exact']);
title('Numerical and exact solution for h');
xlabel('x')
ylabel('h')

% Plotting exact solution for m
subplot(1,2,2)
plot(xx, m_ex(xx, T), 'LineWidth', 1.5, "Color", colors(end, :));
legend_strings = arrayfun(@(x) sprintf('N=%d', x), N_vec, 'UniformOutput', false);
legend([legend_strings, 'exact']);
title('Numerical and exact solution for m');
xlabel('x')
ylabel('m')

% Plotting errors
figure
loglog(err_table(1,:), err_table(2,:), '-o', err_table(1,:), err_table(3,:), '-o', err_table(1,:), err_table(1,:), '-')
legend(["h error", "m error", "dx"])
title('Infinity norm error at final time')
xlabel('dx')
ylabel('error')
