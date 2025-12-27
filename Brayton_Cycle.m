%%  INPUT PARAMETERS
gamma = 1.4;
R = 287;                
cp_const = 1005;        

T1 = 300;               
P1 = 100e3;             
rp = 10;                
T3 = 1400;              

eta_c = 0.85;
eta_t = 0.88;           

%% IDEAL BRAYTON CYCLE

T2 = T1 * rp^((gamma-1)/gamma);
T4 = T3 / rp^((gamma-1)/gamma);

%% NON-IDEAL BRAYTON CYCLE

T2a = T1 * (1 + (T2/T1 - 1)/eta_c);
T4a = T3 * (1 - eta_t*(1 - T4/T3));

%% TEMPERATURE-DEPENDENT cp(T)

cp = @(T) 1000 + 0.1*(T - 300);    

% Enthalpy change using integration
Wc = integral(cp, T1, T2a); % Compressor work
Wt = integral(cp, T4a, T3);% Turbine work

Wnet = Wt - Wc;
Qin  = integral(cp, T2a, T3);

eta_th = Wnet / Qin;
work_ratio = Wnet / Wt;

%% DISPLAING RESULTS

fprintf('\n----- BRAYTON CYCLE RESULTS -----\n');
fprintf('Compressor Work = %.2f J/kg\n', Wc);
fprintf('Turbine Work    = %.2f J/kg\n', Wt);
fprintf('Net Work Output = %.2f J/kg\n', Wnet);
fprintf('Heat Added      = %.2f J/kg\n', Qin);
fprintf('Thermal Efficiency = %.4f\n', eta_th);
fprintf('Work Ratio = %.4f\n', work_ratio);

%% T–s DIAGRAM

s1 = 0;  % reference
s2 = s1 + cp_const*log(T2/T1) - R*log(rp);
s3 = s2 + cp_const*log(T3/T2);
s4 = s3 - cp_const*log(T3/T4);

T_plot = [T1 T2 T3 T4 T1];
s_plot = [s1 s2 s3 s4 s1];

figure;
plot(s_plot, T_plot, '-o', 'LineWidth', 2);
xlabel('Entropy (J/kg·K)');
ylabel('Temperature (K)');
title('T-s Diagram of Brayton Cycle');
grid on;

%% P–v DIAGRAM

v1 = R*T1/P1;
v2 = R*T2/(P1*rp);
v3 = R*T3/(P1*rp);
v4 = R*T4/P1;

P = [P1 P1*rp P1*rp P1];
v = [v1 v2 v3 v4];

figure;
plot(v, P/1e5, '-o', 'LineWidth', 2);
xlabel('Specific Volume (m^3/kg)');
ylabel('Pressure (bar)');
title('P–v Diagram of Brayton Cycle');
grid on;

%% PRESSURE RATIO SWEEP

rp_vec = 2:0.5:20;
Wnet_vec = zeros(size(rp_vec));
eta_vec = zeros(size(rp_vec));

for i = 1:length(rp_vec)
    rp = rp_vec(i);

    T2 = T1 * rp^((gamma-1)/gamma);
    T4 = T3 / rp^((gamma-1)/gamma);

    T2a = T1 * (1 + (T2/T1 - 1)/eta_c);
    T4a = T3 * (1 - eta_t*(1 - T4/T3));

    Wc = integral(cp, T1, T2a);
    Wt = integral(cp, T4a, T3);

    Wnet_vec(i) = Wt - Wc;
    eta_vec(i) = Wnet_vec(i) / integral(cp, T2a, T3);
end

figure;
plot(rp_vec, Wnet_vec, 'LineWidth', 2);
xlabel('Pressure Ratio');
ylabel('Net Work Output (J/kg)');
title('Net Work vs Pressure Ratio');
grid on;

figure;
plot(rp_vec, eta_vec, 'LineWidth', 2);
xlabel('Pressure Ratio');
ylabel('Thermal Efficiency');
title('Efficiency vs Pressure Ratio');
grid on;

