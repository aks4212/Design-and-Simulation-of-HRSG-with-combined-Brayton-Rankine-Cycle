%% IDEAL RANKINE CYCLE

% INPUT PARAMETERS
Pb = 15e6;          
Pc = 10e3;          
T3 = 823;           

% Converting pressures to bar for steam tables
Pb_bar = Pb / 1e5;
Pc_bar = Pc / 1e5;

% STATE 1: CONDENSER EXIT
% Saturated liquid at condenser pressure
h1 = 191.8;    
s1 = 0.649;     
v1 = 0.00101;    

% STATE 2: PUMP EXIT (ISENTROPIC) 
wp = v1 * (Pb - Pc) / 1000;     
h2 = h1 + wp;
s2 = s1;

% STATE 3: BOILER EXIT
T3_C = T3 - 273.15;              
h3 = 3510;
s3 = 6.6;

% STATE 4: TURBINE EXIT (ISENTROPIC) 
s4 = s3;
h4 = 2400;

% PERFORMANCE PARAMETERS 
Wt = h3 - h4;     
Wp = h2 - h1;     
Wnet = Wt - Wp;

Qin = h3 - h2;

eta_th = Wnet / Qin;
BWR = Wp / Wt;

% DISPLAING THE RESULTS 
fprintf('--- IDEAL RANKINE CYCLE ---\n');
fprintf('Turbine Work = %.2f kJ/kg\n', Wt);
fprintf('Pump Work = %.2f kJ/kg\n', Wp);
fprintf('Net Work = %.2f kJ/kg\n', Wnet);
fprintf('Heat Added = %.2f kJ/kg\n', Qin);
fprintf('Thermal Efficiency = %.3f\n', eta_th);
fprintf('Back Work Ratio = %.4f\n', BWR);

%% NON-IDEAL RANKINE CYCLE
eta_t = 0.85;   % Turbine efficiency(let)
eta_p = 0.80;   % Pump efficiency(let)

% Actual pump work
h2_actual = h1 + (h2 - h1)/eta_p;

% Isentropic turbine exit
h4s = h4;

% Actual turbine exit
h4_actual = h3 - eta_t*(h3 - h4s);

% Recalculate performance
Wt_real = h3 - h4_actual;
Wp_real = h2_actual - h1;
Wnet_real = Wt_real - Wp_real;
Qin_real = h3 - h2_actual;

eta_real = Wnet_real / Qin_real;

fprintf('\n--- REAL RANKINE CYCLE ---\n');
fprintf('Net Work = %.2f kJ/kg\n', Wnet_real);
fprintf('Thermal Efficiency = %.3f\n', eta_real);

% plotting
Tsat = 45.8;
T = [Tsat, Tsat, T3_C, Tsat];
s = [s1, s2, s3, s4];
h = [h1, h2, h3, h4];

figure;
plot(s, T, '-o', 'LineWidth', 2);
xlabel('Entropy (kJ/kg·K)');
ylabel('Temperature (°C)');
title('T-s Diagram of Rankine Cycle');
grid on;

figure;
plot(s, h, '-o', 'LineWidth', 2);
xlabel('Entropy (kJ/kg·K)');
ylabel('Enthalpy (kJ/kg)');
title('h-s Diagram of Rankine Cycle');
grid on;

% Boiler Pressure Parametric Study
Pb_range = linspace(5e6, 20e6, 20);
eta_vec = zeros(size(Pb_range));

for i = 1:length(Pb_range)
    Pb = Pb_range(i);

  % taking approximate values from steam tables
    if Pb <= 10e6
        h3 = 3400;
    elseif Pb <= 15e6
        h3 = 3510;
    else
        h3 = 3600;
    end

    % Isentropic turbine outlet (approx)
    h4 = 2400;

    Wt = h3 - h4;
    Wp = v1*(Pb - Pc)/1000;
    Wnet = Wt - Wp;
    Qin = h3 - h2;

    eta_vec(i) = Wnet / Qin;
end

figure;
plot(Pb_range/1e6, eta_vec, 'LineWidth', 2);
xlabel('Boiler Pressure (MPa)');
ylabel('Thermal Efficiency');
title('Effect of Boiler Pressure on Rankine Cycle Efficiency');
grid on;


