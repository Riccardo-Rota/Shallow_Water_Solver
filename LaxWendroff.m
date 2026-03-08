function [xx, t_final, q_final] = LaxWendroff(x_initial, x_final, N, T, bc, h0, m0, f, J, S, CFL, g)
  
    %discretization in space
    xx = linspace(x_initial, x_final, N+1);
    dx = (x_final - x_initial) / N;

    %initial condition
    q = [h0(xx); m0(xx)];
    q_new = nan(2,N+1);

    t=0;
    while t < T
        k = CFL * dx/ max(abs(q(2,:)./q(1,:)) + sqrt(abs(g*q(1,:))));
        for j = 2:N
            q_new(:,j) = q(:,j) - k/(2*dx) * (f(q(:,j+1)) - f(q(:,j-1))) ...
                           + k^2/(2*dx^2) * J((q(:,j+1)+q(:,j))/2) * (f(q(:,j+1)) - f(q(:,j)))...
                           - k^2/(2*dx^2) * J((q(:,j)+q(:,j-1))/2) * (f(q(:,j)) - f(q(:,j-1)));           
    
            if bc == "periodic"
                q_new(:,1) = q(:,1) - k/(2*dx) * (f(q(:,2)) - f(q(:,N))) ...
                           + k^2/(2*dx^2) * J((q(:,2)+q(:,1))/2) * (f(q(:,2)) - f(q(:,1)))...
                           - k^2/(2*dx^2) * J((q(:,1)+q(:,N))/2) * (f(q(:,1)) - f(q(:,N)));
                q_new(:,N+1) = q(:,N+1) - k/(2*dx) * (f(q(:,2)) - f(q(:,N))) ...
                           + k^2/(2*dx^2) * J((q(:,2)+q(:,N+1))/2) * (f(q(:,N+1)) - f(q(:,N*1)))...
                           - k^2/(2*dx^2) * J((q(:,N*1)+q(:,N))/2) * (f(q(:,N+1)) - f(q(:,N)));
        
            elseif bc == "open"
                q_new(:,1) = q(:,1) - k/(2*dx) * (f(q(:,2)) - f(q(:,1))) ...
                           + k^2/(2*dx^2) * J((q(:,2)+q(:,1))/2) * (f(q(:,2)) - f(q(:,1)))...
                           - k^2/(2*dx^2) * J((q(:,1)+q(:,1))/2) * (f(q(:,1)) - f(q(:,1)));
                q_new(:,N+1) = q(:,N+1) - k/(2*dx) * (f(q(:,N+1)) - f(q(:,N))) ...
                           + k^2/(2*dx^2) * J((q(:,N+1)+q(:,N+1))/2) * (f(q(:,N+1)) - f(q(:,N+1)))...
                           - k^2/(2*dx^2) * J((q(:,N+1)+q(:,N))/2) * (f(q(:,N+1)) - f(q(:,N)));
        
            else
                error("Invalid boundary condition");
            end
        end
        t = t + k;
        q = q_new;
    end
    t_final = t;
    q_final = q;

end