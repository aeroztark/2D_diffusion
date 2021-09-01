function[u] = ExplicitMethod(k,CFL,u,dt,dx,dy,jD,iD)

% This function solves 2D diffusion equation using Finite Differencce
% method explicitly. The explicit scheme is forward Euler in time and
% centererd difference in space.

% created by: Sarthak Srivastava, 01-Sep-2021

    % First calculating the number of sub-steps for integration of diffusion equation 
    N_substeps = ceil(dt*k/(CFL*min(dx,dy)^2));   %no of substeps required to solve diffusion equation...
        ... based on Von Neumann Number (dt = N_substeps x dt_sub)
    
    % Main Substepping Loop
    for m = 1:N_substeps
        
        % Substep timestep
        dt_sub = dt/N_substeps;
        
        % x-split for diffusion equation ----
        % Using an explicit scheme: forward Euler in time and centrered difference in space
        u(jD,iD) = u(jD,iD) + k.*(dt_sub/dx^2).*(u(jD,iD+1)-2.*u(jD,iD)+u(jD,iD-1)); % get updated T
        % apply BCs
        u(1,:) = 0;
        u(end,:) = 0;
        u(:,1) = 0;
        u(:,end) = 0;
        
        % z-split for diffusion equation ----
        u(jD,iD) = u(jD,iD) + k.*(dt_sub/dy^2).*(u(jD+1,iD)-2.*u(jD,iD)+u(jD-1,iD)); % get updated T
        % apply BCs
        u(1,:) = 0;
        u(end,:) = 0;
        u(:,1) = 0;
        u(:,end) = 0;
    
    end
end