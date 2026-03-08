%% 1.3

%% Problem setting and exact solution (obtained by Lax-Friedrichs with fine grid)
clc
close
clear
g = 1;
x_initial = 0;
x_final = 2;
T = 0.5;
CFL = 0.25;
f = @(q) [q(2,:) + 0*q(1,:); q(2,:).^2./q(1,:) + 0.5*g*q(1,:).^2]; %Physical flux
J = @(q) [0, 1; -(q(2,:)./q(1,:)).^2 + g*q(1,:), 2*q(2,:)./q(1,:)]; %Jacobian of f
S = @(x,t) [0.*x.*t; 0.*x.*t]; %External source
h0 = @(x) 1 + 0*x; %Initial condition on h
m0 = @(x) -0.5*(x<1); %Initial condition on m
bc="open"; %Setting of the boundary condition (open or periodic)
N_ex = 2^14; %Number of refinements for computing 'exact' solution
N_vec = 2.^[7,8,9,10]; %Number of refinements for the computed numerical solution
colors = lines(length(N_vec)+1); 

%compute the 'exact' solution with LF
[xx_ex, t, q_ex] = LaxFriedrichs(x_initial, x_final, N_ex, T, bc, h0, m0, f, S, CFL, g);
dx= xx_ex(2)-xx_ex(1);

%% Approximate solution with Lax-Friedrichs
figure
%plot of the 'exact' solution
subplot(1,2,1)
plot(xx_ex, q_ex(1,:), 'LineWidth', 1.5, "Color", colors(end, :));
hold on
subplot(1,2,2)
plot(xx_ex, q_ex(2,:), 'LineWidth', 1.5, "Color", colors(end, :));
hold on

err_table = zeros(3, length(N_vec));
for i = 1:length(N_vec)
    N=N_vec(i);
    
    %Compute the numerical solution with LF
    [xx, t, q] = LaxFriedrichs(x_initial, x_final, N, T, bc, h0, m0, f, S, CFL, g);
    dx= xx(2)-xx(1);
    
    %Plot of the numerical solution
    subplot(1,2,1)
    plot(xx, q(1,:), "Color", colors(i, :));
    hold on
    subplot(1,2,2)
    plot(xx, q(2,:), "Color", colors(i, :));
    hold on

    %Computation of the error
    idx_ex= 1: N_ex/N : N_ex+1; %indeces of finer grid corresponding to the coarser grid
    m_err_inf = norm(q(2,:)-q_ex(2,idx_ex), "inf");
    h_err_inf = norm(q(1,:)-q_ex(1,idx_ex), "inf");
    err_table(3,i)=m_err_inf;
    err_table(2,i)=h_err_inf;
    err_table(1,i)=dx;
end

%legend and titles
subplot(1,2,1)
legend_strings = arrayfun(@(x) sprintf('N=%d', x), N_vec, 'UniformOutput', false);
legend(['exact', legend_strings]);
title('Numerical and exact solution for h');
xlabel('x')
ylabel('h')

subplot(1,2,2)
legend_strings = arrayfun(@(x) sprintf('N=%d', x), N_vec, 'UniformOutput', false);
legend(['exact', legend_strings]);
title('Numerical (Lax-f) and exact solution for m');
xlabel('x')
ylabel('m')

%plot of the erorrs
figure
loglog(err_table(1,:), err_table(2,:), '-o', err_table(1,:), err_table(3,:), '-o', err_table(1,:), err_table(1,:), '-')
legend(["h error", "m error", "dx"])
title('Lax-Friedrichs error at final time')
xlabel('dx')
ylabel('error')

%% Approximate solution with Lax-Wendroff
N_vec = 2.^[7,8,9,12]; %For the graphs in the report we change the refinements to these
%values, but it takes some minutes to run, leave them commented for
%a faster plot

%plot of the 'exact' solution
figure
subplot(1,2,1)
plot(xx_ex, q_ex(1,:), 'LineWidth', 1.5, "Color", colors(end, :));
hold on
subplot(1,2,2)
plot(xx_ex, q_ex(2,:), 'LineWidth', 1.5, "Color", colors(end, :));
hold on

err_table = zeros(3, length(N_vec));
for i = 1:length(N_vec)
    N=N_vec(i);
   
    %Compute the numerical solution with LW
    [xx, t, q] = LaxWendroff(x_initial, x_final, N, T, bc, h0, m0, f, J, S, CFL, g);
    dx= xx(2)-xx(1);

    %Plot of the numerical solution
    subplot(1,2,1)
    plot(xx, q(1,:), "Color", colors(i, :));
    hold on
    subplot(1,2,2)
    plot(xx, q(2,:), "Color", colors(i, :));
    hold on

    %Computation of the error
    idx_ex= 1: N_ex/N : N_ex+1; %indeces of finer grid corresponding to the coarser grid
    m_err_inf = norm(q(2,:)-q_ex(2,idx_ex), "inf");
    h_err_inf = norm(q(1,:)-q_ex(1,idx_ex), "inf");
    err_table(3,i)=m_err_inf;
    err_table(2,i)=h_err_inf;
    err_table(1,i)=dx;
end

%legend and titles
subplot(1,2,1)
legend_strings = arrayfun(@(x) sprintf('N=%d', x), N_vec, 'UniformOutput', false);
legend(['exact', legend_strings]);
title('Numerical (Lax-Wendroff) and exact solution for h');
xlabel('x')
ylabel('h')

subplot(1,2,2)
legend_strings = arrayfun(@(x) sprintf('N=%d', x), N_vec, 'UniformOutput', false);
legend(['exact', legend_strings]);
title('Numerical (Lax-Wendroff) and exact solution for m');
xlabel('x')
ylabel('m')

%plot of the erorrs
figure
loglog(err_table(1,:), err_table(2,:), '-o', err_table(1,:), err_table(3,:), '-o', err_table(1,:), err_table(1,:), '-')
legend(["h error", "m error", "dx"])
title('Lax-Wendroff error at final time')
xlabel('dx')
ylabel('error')