% Material data relevant to the structure

TI = 1; % Titanium
GL = 0; % Glass

D = containers.Map({TI,GL},{eye(2)*17, eye(2)*0.8}); % Thermal conductivity [W/(m K)] 
alpha_c = 100; % Heat transfer coefficient [W/(m^2K)]
T_0 = 20; % Temperature at which structure is stress free [C]

rho = containers.Map({TI,GL},{4620,3860}); % Density [kg/m^3]
E = containers.Map({TI,GL},{110, 67}); % Young's modulus [GPa]
alpha = containers.Map({TI,GL},{9.4e-6,7e-6}); % Expansion coefficient [1/K]
c_p = containers.Map({TI,GL},{523, 670}); % Specific heat [J/(kg K)]
Poisson = containers.Map({TI,GL},{0.34, 0.2}); % Poission's ratio [-]

thickness = 0.01; % Thickness of disc structure [m]