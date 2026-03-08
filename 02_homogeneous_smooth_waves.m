%% 1.2

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
S = @(x,t) [0.*x.*t; 0.*x.*t]; % source term
% Two possible initial conditions
h0_1 = @(x) 1 - 0.1*sin(pi*x);
m0_1 = @(x) 0.*x;
h0_2 = @(x) 1 - 0.2*sin(2*pi*x);
m0_2 = @(x) 0.5+ 0.*x;

% CHOOSE HERE which initial condition consider!
h0 = h0_2;
m0 = m0_2;

% Set typer of boundary conditions
bc="periodic";
%Very fine mesh for "exact" solution
N_ex = 2^14;
% Mesh for numerical solutions
N_vec = 2.^[9,10,11,12];
% Initilize error vector and colors for plots
err_table = zeros(3, length(N_vec));
colors = lines(length(N_vec)+1);

% Create exact solution
[xx_ex, t, q_ex] = LaxFriedrichs(x_initial, x_final, N_ex, T, bc, h0, m0, f, S, CFL, g);
dx= xx_ex(2)-xx_ex(1);

figure
subplot(1,2,1)
plot(xx_ex, q_ex(1,:), 'LineWidth', 1.5, "Color", colors(end, :));
hold on
subplot(1,2,2)
plot(xx_ex, q_ex(2,:), 'LineWidth', 1.5, "Color", colors(end, :));
hold on

% Compute numerical solution for each grid
for i = 1:length(N_vec)
    N=N_vec(i);
   
    [xx, t, q] = LaxFriedrichs(x_initial, x_final, N, T, bc, h0, m0, f, S, CFL, g);
    dx= xx(2)-xx(1);

    subplot(1,2,1)
    plot(xx, q(1,:), "Color", colors(i, :));
    hold on
    subplot(1,2,2)
    plot(xx, q(2,:), "Color", colors(i, :));
    hold on
    % Compute and save errors (peak the right points of exact solution grid using
    % idx_ex)
    idx_ex= 1: N_ex/N : N_ex+1;
    m_err_inf = norm(q(2,:)-q_ex(2,idx_ex), "inf");
    h_err_inf = norm(q(1,:)-q_ex(1,idx_ex), "inf");
    err_table(3,i)=m_err_inf;
    err_table(2,i)=h_err_inf;
    err_table(1,i)=dx;
end

% Plot all the results
subplot(1,2,1)
legend_strings = arrayfun(@(x) sprintf('N=%d', x), N_vec, 'UniformOutput', false);
legend(['exact', legend_strings]);
title('Numerical and exact solution for h');
xlabel('x')
ylabel('h')

subplot(1,2,2)
legend_strings = arrayfun(@(x) sprintf('N=%d', x), N_vec, 'UniformOutput', false);
legend(['exact', legend_strings]);
title('Numerical and exact solution for m');
xlabel('x')
ylabel('m')

figure
loglog(err_table(1,:), err_table(2,:), '-o', err_table(1,:), err_table(3,:), '-o', err_table(1,:), err_table(1,:), '-')
legend(["h error", "m error", "dx"])
title('Infinity norm error at final time')
xlabel('dx')
ylabel('error')