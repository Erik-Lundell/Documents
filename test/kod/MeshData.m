% Extract CALFEM notation from pdetool-mesh

% Element data
% Row #i of matrix refers to element #i

enod=t(1:3,:)'; % Nodes of elements
emat = (t(4,:)==TI)';  %  Material type of element
nelm = size(enod,1); % Number of elements

% Node data

coord = p'/100; % Coordinates of each node
nnod = size(coord,1); % Number of nodes
dof = (1:nnod)'; % Degrees of freedom
dof_S = [(1:nnod)',(nnod+1:2*nnod)']; % Give each dof a number
edof_S = zeros(nelm, 7);   
edof = zeros(nelm, 4); %Element degrees of freedom

for ie = 1:nelm
    edof_S(ie,:) = [ie dof_S(enod(ie,1),:), dof_S(enod(ie,2),:),dof_S(enod(ie,3),:)];
    edof(ie,:) = [ie,enod(ie,:)];
end

% Define boundary conditions of element edges
% An edge of an element is defined by a node pair, [node #1; node #2]

er = e([1 2 5],:); % [edge; boundary segment]
conv_segments = [2 15 14 23 31 1]; % Boundary segments with convection
edges_conv = []; % Edges with convection

fixed_segments = [21 22 1 16 17 18 19 20]; % Fixed boundary segments
edges_fixed = []; % Fixed edges

for i = 1:size(er,2)
    if ismember(er(3,i),conv_segments)
        edges_conv = [edges_conv [er(1:2,i);er(3,i)==1]]; % [edge; convection type]
            % convection type = 0 -> T_inf outside
            % convection type = 1 -> T_c outside
    end
    if ismember(er(3,i),fixed_segments)
        edges_fixed = [edges_fixed [er(1:2,i);er(3,i)==1]]; % [edge; fixed type]
            %fixed type = 0 -> u_y = 0
            %fixed type = 1 -> u_x = 0
    end
end

