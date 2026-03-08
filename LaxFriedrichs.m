function [xx, t_final, q_final] = LaxFriedrichs(x_initial, x_final, N, T, bc, h0, m0, f, S, CFL, g)
  
    %discretization in space
    xx = linspace(x_initial, x_final, N+1);
    dx = (x_final - x_initial) / N;

    %initial condition
    q = [h0(xx); m0(xx)];
    q_new = nan(2,N+1);

    t=0;
    while t < T
        k = CFL * dx/ max(abs(q(2,:)./q(1,:)) + sqrt(abs(g*q(1,:))));
        alpha=dx/k;
        q_new(:,2:N) = q(:,2:N) - k/(2*dx) * (f(q(:,3:N+1)) - f(q(:,1:N-1))) + ...
                       alpha*k/(2*dx)*(q(:,3:N+1) - 2*q(:,2:N) + q(:,1:N-1)) + k*S(xx(2:N), t);
        if bc == "periodic"
            q_new(:,1) = q(:,1) - k/(2*dx) * (f(q(:,2)) - f(q(:,N))) + ...
                       alpha*k/(2*dx)*(q(:,2) - 2*q(:,1) + q(:,N)) + k*S(xx(1), t);
            q_new(:,N+1) = q(:,N+1) - k/(2*dx) * (f(q(:,2)) - f(q(:,N))) + ...
                       alpha*k/(2*dx)*(q(:,2) - 2*q(:,N+1) + q(:,N)) + k*S(xx(N+1), t);
    
        elseif bc == "open"
            q_new(:,1) = q(:,1) - k/(2*dx) * (f(q(:,2)) - f(q(:,1))) + ...
                     alpha*k/(2*dx)*(q(:,2) - 2*q(:,1) + q(:,1)) + k*S(xx(1), t);
            q_new(:,N+1) = q(:,N+1) - k/(2*dx) * (f(q(:,N+1)) - f(q(:,N))) + ...
                     alpha*k/(2*dx)*(q(:,N+1) - 2*q(:,N+1) + q(:,N)) + k*S(xx(N+1), t);
    
        else
            error("Invalid boundary condition");
        end
        t = t + k;
        q = q_new;
    end
    t_final = t;
    q_final = q;

end