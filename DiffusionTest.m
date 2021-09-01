% This script compares the Explicit and Implicit numerical methods against the exact
% solution of a 2-D sine hill diffusion problem

% created by: Sarthak Srivastava, 01-Sep-2021

clear
clc

% Setting up the domain
k = 100;
Lx = 100;
Ly = 100;
dx = 1;
dy = 1;
param = -(k*pi^2)/(Lx^2 + Ly^2); % for ease of computation
x = 1:dx:Lx;
y = 1:dy:Ly;

% time parameters
t = 0;
dt = 0.8;
CFL = 0.5; % for Explicit method
T_end = 50;
t_array = 0:dt:T_end;

% array init
max_err = zeros(size(t_array));
rms_err = zeros(size(t_array));
u_exact = zeros(length(x),length(y));
u_explicit = zeros(length(x),length(y));
u_implicit = zeros(length(x),length(y));

% Initial Conditions
for i = 1:Lx
    for j = 1:Ly
        u_exact(i,j) = 10.*sin((pi/Lx).*x(i)).*sin((pi/Ly).*y(j)) ;
        u_implicit(i,j) = 10.*sin((pi/Lx).*x(i)).*sin((pi/Ly).*y(j)) ;
        u_explicit(i,j) = 10.*sin((pi/Lx).*x(i)).*sin((pi/Ly).*y(j)) ;
    end
end

figure
set(gcf,'Position',[100,500,1800,400])

%% time loop
for ind = 1:length(t_array)
    
% Set boundary conditions
    u_exact(1,:) = 0;
    u_exact(end,:) = 0;
    u_exact(:,1) = 0;
    u_exact(:,end) = 0;
      
% Exact solution
    for i = 1:Lx
        for j = 1:Ly
            u_exact(i,j) = 10.*exp(param*t).*sin((pi/Lx).*x(i)).*sin((pi/Ly).*y(j)) ;
        end
    end

% Solution by Explicit method
    u_explicit = ExplicitMethod(k,CFL,u_explicit,dt,dx,dy,2:length(y)-1,2:length(x)-1);
    Explicitdiff = (u_exact - u_explicit)./10;
    Explicitmax_err(ind) = max(max(abs(Explicitdiff)));
    Explicitrms_err(ind) = rms(reshape(Explicitdiff,[length(x)*length(y),1]));
    
 % Solution by Implicit method
    u_implicit = ImplicitMethod(k,u_implicit,dx,length(x),length(y));
    Implicitdiff = (u_exact - u_implicit)./10;
    Implicitmax_err(ind) = max(max(abs(Implicitdiff)));
    Implicitrms_err(ind) = rms(reshape(Implicitdiff,[length(x)*length(y),1]));
    
 % Visualizing the solution   
    subplot(1,3,1)
    surf(u_exact)
    zlim([0 10])
    title('Exact')
    colorbar
    subplot(1,3,2)
    surf(u_explicit)
    zlim([0 10])
    title('Explicit')
    colorbar
    subplot(1,3,3)
    surf(u_implicit)
    zlim([0 10])
    title('Implicit')
    colorbar
    pause(0.01)
    ind = ind+1;
    t = t+dt;
end

% Scaled error plot
figure
subplot(1,2,1)
plot(t_array,Implicitmax_err)
hold on
plot(t_array,Explicitmax_err)
legend('Implicit','Explicit')
xlabel('timestep')
ylabel('Scaled max error')
subplot(1,2,2)
plot(t_array,Implicitrms_err)
hold on
plot(t_array,Explicitrms_err)
legend('Implicit','Explicit')
xlabel('timestep')
ylabel('Scaled RMS error')
