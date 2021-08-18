% Transient heat problem

delta_t = 100; % Time step length
T_outside = [-96 20]; % Initial boundary temperatures [T_inf T_c]
a = a0; % Initial nodal temperatures, generated in a)

A = zeros(nnod); % Damping matrix
F = zeros(nnod,1); % Force vector

% Generate F with initial boundary temperatures

fce_const = [1; 1]*alpha_c/2;
for ib = 1:length(edges_conv)
    indx = edges_conv(1:2,ib); % Edge with convection
    ex = coord(indx,1); % x-coordinates of nodes
    ey = coord(indx,2); % y-coordinates of nodes
    l = sqrt((ex(1)-ex(2))^2+(ey(1)-ey(2))^2);  % Edge length
    
    convectionTemp = edges_conv(3,ib)+1; % Nearby outside temperature 
    fce = fce_const*l*T_outside(convectionTemp); % Force vector integral
    F(indx) = F(indx) + fce; % Insert into F
end

% Damping matrix A

for ie = 1:nelm
    ex = coord(enod(ie,:),1); % x-coordinates of nodes
    ey = coord(enod(ie,:),2); % y-coordinates of nodes
   
    material = emat(ie); % Material type of element
    Ae = plantml(ex',ey',rho(material)*c_p(material)); % Element damping matrix
    
    indx = edof(ie,2:end);  % Position in A
    A(indx, indx) = A(indx,indx)+Ae; % Insert into A
end

% Implicit Euler time step function

time_step = @(a) (A + delta_t*K)\(F*delta_t+A*a);

% Temperature development over time

for i=1:200 % 200 time steps
    a = time_step(a);
end

[ex, ey] = coordxtr(edof, coord, dof, 3); % Extract nodal coordinate data
ed = extract(edof, a); % Extract element temperatures

%Plot

patch(ex',ey',ed','EdgeColor','none');
hold on
patch(ex',-ey',ed','EdgeColor','none');

caxis([-100 50]);
axis([-0.1 1.2 -0.5 0.5]/100);
title("t = " + (i*delta_t) + " s");
colormap default;
colorbar;
xlabel('x-position [m]');
ylabel('y-postition [m]');

disp("MAX TEMP: " + max(max(ed)));
disp("MIN TEMP: " + min(min(ed)));