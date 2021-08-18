% Stationary temporary problem

T_outside = [-96 20]; % Nighttime boundary temperatures [T_inf T_c]

K = zeros(nnod); % Global stiffness matrix
F = zeros(nnod,1); % Force vector

% Basic stiffness matrix

for ie = 1:nelm
    ex = coord(enod(ie,:),1)'; % x-coordinates of element
    ey = coord(enod(ie,:),2)'; % y-coordinates of element

    Ke = flw2te(ex,ey,thickness,D(emat(ie))); % Basic element stiffness matrix
    
    indx = edof(ie,2:end);  % Position in K
    K(indx,indx) = K(indx,indx)+Ke;  % Insert into K
end

% Add contribution to K and F from convection

fce_const = [1; 1]*alpha_c/2; % Common force vector integral constant 
Kce_const = [2 1; 1 2]*alpha_c/6; % Common stiffness matrix integral constant
    
for ib = 1:length(edges_conv)
    indx = edges_conv(1:2,ib); % Edge with convection
    ex = coord(indx,1); % x-coordinates of nodes
    ey = coord(indx,2); % y-coordinates of nodes
    l = sqrt((ex(1)-ex(2))^2+(ey(1)-ey(2))^2);  % Edge length
    
    convectionTemp = edges_conv(3,ib)+1; % Nearby outside temperature 
    fce = fce_const*l*T_outside(convectionTemp); % Force vector integral
    Kce = Kce_const*l; % Stiffness matrix integral

    K(indx,indx) = K(indx,indx)+Kce;  % Insert into K
    F(indx) = F(indx) + fce; % Insert into F
end

% Solve system

a0 = solveq(K,F); % Nodal temperatures
[ex, ey] = coordxtr(edof, coord, dof, 3); % Extract nodal coordinate data
ed = extract(edof, a0); % Extract element temperatures

% Plot

clf
patch(ex',ey',ed','EdgeColor','none');
hold on
patch(ex',-ey',ed','EdgeColor','none');

caxis([-100 50]);
axis([-0.1 1.2 -0.5 0.5]/100);
title("Temperature distribution, T_{\infty} = " + (T_outside(1)) + " [C]");
%colormap(hot);
colorbar;
xlabel('x-position [m]')
ylabel('y-postition [m]');

disp("MAX TEMP: " + max(max(ed)));

% Max T vid T_inf = -96: 19,9769
% Max T vid T_inf = 40: 40,0093