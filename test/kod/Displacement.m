lens_subdomain = 3;
lens_displacement = 0; % Total lens displacement

for ie = 1:nelm
    if(t(4,ie)==lens_subdomain) % if element is within lens subdomain...
        ex = coord(enod(ie,:),1)'; % x-coordinates of nodes
        ey = coord(enod(ie,:),2)'; % y-coordinates of nodes
        
        T = NtN(ex, ey,thickness); % Compute int(N^T*N)t dA
        a_x = a_S(edof_S(ie,2:4)); % Nodal x-component displacement
        a_y = a_S(edof_S(ie,5:7)); % Nodal y-component displacement
        a = [a_x; a_y]; % Nodal displacement
        
        lens_displacement = lens_displacement + a'*T*a; % Compute lens displacement and add
    end
end

disp("TOTAL LENS DISPLACEMENT: " + lens_displacement);

% Plot

mag = 10; % Displacement magnification
exd = ex_S + mag*ed_S(:,1:2:end); % Nodal x-coordinate data
eyd = ey_S + mag*ed_S(:,2:2:end); % Nodal y-coordinate data

figure();
patch(ex_S',ey_S',[0 0 0],"EdgeColor","none","FaceAlpha",0.3);
hold on
patch(ex_S',-ey_S',[0 0 0],"EdgeColor","none","FaceAlpha",0.3);
patch(exd',eyd',[0 0 0],"FaceAlpha",0.3);
patch(exd',-eyd',[0 0 0],"FaceAlpha",0.3);
axis equal

xlabel('x-position [m]');
ylabel('y-postition [m]');
title("Node displacements (magnified x10)");