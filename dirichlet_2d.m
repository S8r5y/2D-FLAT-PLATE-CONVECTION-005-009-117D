% Parameters
L = 1;          
W = 0.1;        
U_infty = 0.1;  
T_plate = 100;                                                                                                      
T_infty = 20;   
k = 0.1;        
h = 0.01;       
Q_gen = 1000;   

% Grid parameters
Nx = 20;        % Number of grid points in x-direction
Ny = 10;        % Number of grid points in y-direction
dx = L / Nx;    % Grid spacing in x-direction
dy = W / Ny;    % Grid spacing in y-direction

% Time parameters
dt = 0.1;       % Time step
t_end = 100;    % End time
Nt = round(t_end / dt); % Number of time steps

% Thermal properties
alpha = k / (rho * cp); % Thermal diffusivity

% Initialize temperature array
T = T_infty * ones(Nx+1, Ny+1);

% Main time loop
for n = 1:Nt
    % Copy temperature array for the next time step
    T_new = T;
    
    % Compute temperature distribution using explicit finite difference method
    for i = 2:Nx
        for j = 2:Ny
            % Compute Laplacian using central differences
            laplacian_T = (T(i+1,j) - 2*T(i,j) + T(i-1,j)) / dx^2 + (T(i,j+1) - 2*T(i,j) + T(i,j-1)) / dy^2;
            
            % Compute advection term
            adv_T = -U_infty * ((T(i,j) - T(i-1,j)) / dx);
            
            % Update temperature
            T_new(i, j) = T(i, j) + dt * (alpha * laplacian_T + adv_T + Q_gen / (rho * cp));
        end
    end
    
    % Update temperature array for the next time step
    T = T_new;
    
    % Apply Dirichlet boundary conditions
    T(1,:) = T_plate;
    T(Nx+1,:) = T_infty;
    T(:,1) = T_infty;
    T(:,Ny+1) = T_infty;
end

% Plot temperature distribution
[X, Y] = meshgrid(0:dx:L, 0:dy:W);
contourf(X, Y, T', 'LineColor', 'none');
colorbar;
xlabel('X');
ylabel('Y');
title('Temperature Distribution (2D Convection with Heat Generation)');
