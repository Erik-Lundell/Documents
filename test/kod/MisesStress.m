% Mechanical problem 

% D matrix for isotropic material with plane strain

calcD = @(E, v) E/(1+v)*(1-2*v)*[(1-v) v 0; v (1-v) 0; 0 0 (1-2*v)/2];
D_el = containers.Map({TI,GL},{calcD(E(TI),Poisson(TI)),calcD(E(TI),Poisson(GL))});

% Define constant part of D*epsilon_0

calcConstEps_0 = @(mat) alpha(mat)*E(mat)/(1-2*Poisson(mat))*[1;1;0];
const_Deps0 = containers.Map({TI,GL},{calcConstEps_0(TI),calcConstEps_0(GL)});

K = zeros(nnod*2); % Global stiffness matrix
F = zeros(nnod*2,1); % Force vector

for ie = 1:nelm
    ex = coord(enod(ie,:),1)'; % x-coordinates of nodes
    ey = coord(enod(ie,:),2)'; % y-coordinates of nodes
    material = emat(ie); % Material type of element
    
    Ke = plante(ex, ey, [2 thickness], D_el(material)); % Element stiffness matrix
    
    dT = mean(ed(ie,:))-T_0; % Mean temperature in element
    es = const_Deps0(material)*dT; % Element D*epsilon_0
    f_0e = plantf(ex, ey, [2 thickness], es'); % Element f_0
    % (f_be = 0)
    
    indx = edof_S(ie,2:end);  % Position in matrix
    K(indx, indx) = K(indx,indx)+Ke; % Insert into K
    F(indx) = F(indx) + f_0e; % Insert into F
end

% Calculate bc

already_added = zeros(nnod*2,1); % Memory vector
bc = []; % Boundrary condition vector

for ib = 1:length(edges_fixed)
    edge = edges_fixed(:,ib); % [edge; fixed type]
    x_or_y = 2-edge(3); % Fixed type
    
    %Check first node
    
    node_id = dof_S(edge(1),x_or_y); % Index of fixed node component
    
    if(already_added(node_id)==0) % if it hasn't been added...
        bc = [bc; node_id, 0]; % Add to bc
        already_added(node_id) = 1; % Update memory vector
    end
    
    % Check second node
    
    node_id = dof_S(edge(2),x_or_y); % Index of fixed node component
    
    if(already_added(node_id)==0) % if it hasn't been added...
        bc = [bc; node_id, 0]; % Add to bc
        already_added(node_id) = 1; % Update memory vector
    end
end

% Solve system

a_S = solveq(K,F, bc); % Nodal component displacements
[ex_S, ey_S] = coordxtr(edof_S, coord, dof_S, 3); % Extract nodal coordinate data
ed_S = extract(edof_S, a_S); % Extract element displacements

% Calculate von Mises stress per element

Seff_el = zeros(nelm,1); % Element von Mises stress

for ie = 1: nelm
    ex = coord(enod(ie,:),1)'; % x-coordinates of nodes
    ey = coord(enod(ie,:),2)'; % y-coordinates of nodes
    material = emat(ie); % Material type of element
    dT = mean(ed(ie,:))-T_0; % Mean temperature in element
    
    a_index = [dof_S(enod(ie,:), 1) ; dof_S(enod(ie,:), 2)]; % Nodal component indices in a_S
    
    % [sigma_xx sigma_yy sigma_xy]
    sigma1 = plants(ex, ey, [2 thickness], D_el(material),a_S(a_index)'); % No initial strain
    sigma = sigma1 - (const_Deps0(material)*dT)'; % Add contribution from initial strain
    
    % sigma_zz
    sigma_zz1 = Poisson(material)*(sigma(1) + sigma(2)); % No initial strain
    sigma_zz = sigma_zz1 - alpha(material)*E(material)*dT/(1-2*Poisson(material)); % Add contribution from initial strain
    
    vonMisesSquared = sigma*sigma' + sigma_zz^2 - sigma(1)*sigma(2)-sigma(1)*sigma_zz-sigma(2)*sigma_zz+2*sigma(3)^2;
    Seff_el(ie) = sqrt(vonMisesSquared);
end

% Calculate nodal von Mieses stress as mean of connected elements

Seff_nod = zeros(nnod,1); % Nodal von Mises stress

for i=1:nnod
    [c0, c1] = find(edof(:,2:4)==i); % Row indices of connected elements
    Seff_nod(i,1) = sum(Seff_el(c0))/size(c0,1); % Mean von Mises stress
end

eM = extract(edof, Seff_nod); % Extract nodal von Mises stress

% Plot

patch(ex_S',ey_S',eM','EdgeColor','none');
hold on
patch(ex_S',-ey_S',eM','EdgeColor','none');

axis([-0.1 1.2 -0.5 0.5]/100);
title("Von Mises effective stress field.");
colormap default;
colorbar;
xlabel('x-position [m]');
ylabel('y-postition [m]');