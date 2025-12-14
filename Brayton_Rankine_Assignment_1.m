%% Problem 1: Internal Energy Estimation With Non-Constant Heat Capacity

% C(T) = 700 + 0.35T + -2x10^-4T^2 J/kg-K
% T1 = 320K
% T2 = 820K

T1 = 320;
T2 = 820;

C_T = @(T) 700 + 0.35.*T - 2e-4.*T.^2;

% trapz method(Numwrical Integration)
N_trapz = 1000;
T_trapz = linspace(T1, T2, N_trapz+1);
C_values = C_T(T_trapz);
delta_U_trapz = trapz(T_trapz,C_values);

% integral method
delta_U_integral = integral(C_T, T1, T2);

% Symbolic Integration
syms T
C_sym = 700 + 0.35*T - 2e-4*T^2;
delta_U_sym = int(C_sym, T, T1, T2);
delta_U_analytical = double(delta_U_sym);

fprintf('P1:\n');
fprintf(' Numerical(trapz): %.4f J/kg\n', delta_U_trapz);
fprintf(' Numerical(integral): %.4f J/kg\n', delta_U_integral);
fprintf(' Analytical(Symbolic): %.4f J/kg\n', delta_U_analytical);

% Determining the minimum Grid Resolution :-
target_err_percentage = 0.1;
targe_err_val = target_err_percentage/100 * delta_U_analytical;

N_min = 1;
curr_err = inf;

while curr_err > targe_err_val
    N_min = N_min + 1;
    T_test = linspace(T1, T2, N_min+1);
    delta_U_test = trapz(T_test, C_T(T_test));
    curr_err = abs(delta_U_test - delta_U_analytical);
end

fprintf(' Minimum grid resolution(N) for error < 0.1%%: %d intervals\n', N_min);
fprintf('\n');

%% Problem 2: Energy Balance for a Closed System with Time-Varying Heating

% u(T) = 450 + 1.1T + 0.0012T^2 (kJ/kg)
% du/dt = q_ext(t) + r(T)
% with q_ext(T) = 5000e^-0.002t kJ/s and r(T) = 1500(1 - e^-0.01T)

T0 = 300;
u0 = 450 + 1.1*T0 + 0.0012*T0^2;

t_span = [0, 4000];

% writing T as a function of u
T_of_u = @(u) (-1.1 + sqrt(1.1^2 - 4*0.0012*(450-u))) / (2*0.0012);

dudt = @(t, u) 5000*exp(-0.002*t) + 1500*(1 - exp(-0.01*T_of_u(u)));

[t, u] = ode45(dudt, t_span, u0);

%converting u(t) to T(t) 
T = T_of_u(u);

q_ext = 5000*exp(-0.002*t);
r_T = 1500*(1 - exp(-0.01*T));

% Finding the time where r(T) > q_ext(t) for the first time
idx_crossover = find(r_T > q_ext, 1, 'first');

if ~isempty(idx_crossover)
    t_crossover = t(idx_crossover);
else
    t_crossover = NaN;
end

fprintf('P2:\n');
fprintf('  Internal energy profile: u(t) calculated. (Plot)\n');
fprintf('  Temperature profile: T(t) calculated. (Plot)\n');
if isnan(t_crossover)
    fprintf('  Reaction heat never surpasses external heating in the interval.\n');
else
    fprintf('  Reaction heat surpasses external heating at approximately t = %.2f s\n', t_crossover);
end

% Plotting Results
figure(1);
subplot(2,1,1);
plot(t, u);
title('Problem 2: Internal Energy Profile u(t)');
xlabel('Time (s)');
ylabel('Internal Energy u (kJ/kg)');
grid on;

subplot(2,1,2);
plot(t, T);
title('Problem 2: Temperature Profile T(t)');
xlabel('Time (s)');
ylabel('Temperature T (K)');
grid on;
fprintf('\n');

%% Problem 3: Entropy in a Compression Process With Real-Gas Effects

PA = 1; 
TA = 300;
PB = 20; 
n = 1.28;
cp = 1.05;
R = 0.287; 

% Determining T(P) along the polytropic path.
TB = TA * (PB/PA)^((n-1)/n);
fprintf('P3:\n');
fprintf('  Final temperature T_B along polytropic path: %.2f K\n', TB);

% Setting up the integration path.
P_grid = linspace(PA, PB, 1000);
T_grid = TA * (P_grid/PA).^((n-1)/n);

Z = @(T, P) 1 + 0.0008.*P./T;

dTdP = @(P) ((n-1)/n) * (TA/PA^((n-1)/n)) * P.^((n-1)/n - 1); 

ds_dP = (1./P_grid) .* (cp * ((n-1)/n) - R * Z(T_grid, P_grid));

% Numerical Entropy Change (Real Gas)
delta_s_real = trapz(P_grid, ds_dP);
fprintf('  Numerical Entropy Change (Real Gas Delta s): %.4f kJ/kg.K\n', delta_s_real);

% Ideal-Gas Entropy Change (Z=1)
ds_dP_ideal = (1./P_grid) .* (cp * ((n-1)/n) - R);
delta_s_ideal = trapz(P_grid, ds_dP_ideal);
fprintf('  Ideal-Gas Entropy Change (Delta s_(ideal)): %.4f kJ/kg.K\n', delta_s_ideal);

% percent deviation
percent_deviation = abs((delta_s_real - delta_s_ideal) / delta_s_ideal) * 100;
fprintf('  Percent Deviation: %.2f%%\n', percent_deviation);
fprintf('  Real gas effects, primarily captured by the Z factor (which introduces P/T dependence),\n');
fprintf('  cause a deviation from the ideal gas entropy change. In this case, Z > 1, which\n');
fprintf('  generally leads to a smaller (more negative) entropy change for the compression process.\n');
fprintf('\n');

%% Problem 4: Internal Energy and Entropy Tracking in an Open System

T1_P4 = 310; T2_P4 = 670; Tb = T2_P4;
h = @(T) 300 + 2.5*T + 0.0007*T.^2;
s = @(T) 2.0*log(T) + 0.001*T;

delta_h = h(T2_P4) - h(T1_P4);
delta_s = s(T2_P4) - s(T1_P4);

Q_grid = linspace(20, 100, 50);
m_grid = linspace(min(Q_grid)/delta_h * 0.5, max(Q_grid)/delta_h * 1.5, 50);
[Q_mesh, m_mesh] = meshgrid(Q_grid, m_grid);
S_gen_mesh = m_mesh * delta_s - Q_mesh / Tb;
is_feasible = S_gen_mesh >= 0;

fprintf('P4: Delta_h: %.2f kJ/kg. Feasible region plotted.\n', delta_h);

figure(2);
contourf(Q_mesh, m_mesh, is_feasible, [0.5 0.5]);
colormap([1 1 1; 0.6 0.9 0.6]); 
hold on;
plot(Q_grid, Q_grid / delta_h, 'k--', 'LineWidth', 2); 
title('P4: Feasible Operating Region (delta-S-gen > 0)');
xlabel('Heat Input Q-dot (kW)');
ylabel('Mass Flow Rate m-dot (kg/s)');
legend('Feasible Region', 'Q-dot = m-dot/delta-h line', 'Location', 'NorthWest');
grid on; hold off;
fprintf('\n');

%% Problem 5: ExergyAnalysisWithTemperature-Dependent Properties

T1_P5 = 350; T2_P5 = 900; T0 = 298;
cp_P5 = @(T) 1200 + 0.4.*T - 1.2e-4.*T.^2; 

delta_s_P5 = integral(@(T) cp_P5(T) ./ T, T1_P5, T2_P5);
delta_h_P5 = integral(cp_P5, T1_P5, T2_P5);

I_values = [0.02, 0.10];
X_dest_values = I_values * delta_h_P5;

I_plot = linspace(0, 0.20, 100);
X_dest_plot = I_plot * delta_h_P5;

fprintf('P5: Delta-s: %.4f J/kg.K. X-dot-dest=%.2f J/kg.\n', delta_s_P5, X_dest_values(1));

figure(3);
plot(I_plot * 100, X_dest_plot);
title('P5: Exergy Destruction vs. Irreversibility Level');
xlabel('Irreversibility Level I (%)');
ylabel('Exergy Destruction X-dot-dest (J/kg)');
grid on;

fprintf('\n');

%% Problem 6: A Fully Dynamic Energy–Entropy Cycle Simulation

t_span_P6 = [0, 5000];
t_grid_P6 = linspace(t_span_P6(1), t_span_P6(2), 5000);

Th = @(t) 900 - 300*exp(-0.0008*t);
Tc = @(t) 300 + 40*sin(0.002*t);
Q_in = @(t) 20000*(1 + 0.3*sin(0.003*t));
eta = @(t) 1 - Tc(t) ./ Th(t);
P = @(t) eta(t) .* Q_in(t);

P_values = P(t_grid_P6);
W_t = cumtrapz(t_grid_P6, P_values);
W_total = trapz(t_grid_P6, P_values);
S_gen_dot = Q_in(t_grid_P6) .* (1./Tc(t_grid_P6) - 1./Th(t_grid_P6));

fprintf('P6: Total Work: %.2f kJ. Max eta: %.4f\n', W_total, max(eta(t_grid_P6)));

figure(4);
subplot(3,1,1); plot(t_grid_P6, eta(t_grid_P6) * 100); title('P6: Efficiency eta(t)'); xlabel('Time (s)'); ylabel('Efficiency (%)'); grid on;
subplot(3,1,2); plot(t_grid_P6, W_t); title('P6: Cumulative Work Output W(t)'); xlabel('Time (s)'); ylabel('Work W (kJ)'); grid on;
subplot(3,1,3); plot(t_grid_P6, S_gen_dot); title('P6: Entropy Generation Rate S-dot-gen(t)'); xlabel('Time (s)'); ylabel('S-dot-gen (kW/K)'); grid on;

fprintf('\n');

%% Problem 7: Polytropic Piston With Temperature-Dependent Internal Energy

PA_P7 = 1; TA_P7 = 300; PB_P7 = 10; m_P7 = 1.25; R_P7 = 0.287;

TB_P7 = TA_P7 * (PB_P7/PA_P7)^((m_P7-1)/m_P7);
u = @(T) 500 + 0.8*T + 1.5e-3*T.^2;
delta_u_P7 = u(TB_P7) - u(TA_P7);

W_analytical = R_P7 * (TB_P7 - TA_P7) / (1 - m_P7);
Q = delta_u_P7 + W_analytical;

VA = R_P7 * TA_P7 / PA_P7;
V_of_P = @(P) VA * (P/PA_P7).^(-1/m_P7);
V_grid_P7_num = linspace(V_of_P(PA_P7), V_of_P(PB_P7), 1000);
P_grid_P7_num = PA_P7 * (V_grid_P7_num / VA).^(-m_P7);
W_numerical = -trapz(V_grid_P7_num * 100, P_grid_P7_num * 10);

percent_dev_W = abs((W_numerical - W_analytical) / W_analytical) * 100;
fprintf('P7: TB: %.2f K. Delta-U: %.2f, W: %.2f, Q: %.2f kJ/kg. W Dev: %.2f%%\n',...
    TB_P7, delta_u_P7, W_analytical, Q, percent_dev_W);

fprintf('\n');

%% Problem 8: Optimal Heating Profile Minimizing Entropy Generation

c = 2.5; T1_P8 = 300; Q_total = 5e5; Tb_P8 = 300;
N_q = 50; tf = 2000; dt = tf / N_q; q_initial = Q_total / tf;

objective_entropy = @(q) objective_entropy_P8(q, T1_P8, Tb_P8, c, dt);
Aeq = ones(1, N_q); beq = Q_total / dt; lb = zeros(N_q, 1);
options = optimoptions('fmincon', 'Display', 'off');

[q_optimal, S_gen_min_rate_sum] = fmincon(objective_entropy, q_initial * ones(N_q, 1), [], [], Aeq, beq, lb, [], [], options);
S_gen_optimal = S_gen_min_rate_sum * dt;

q_uniform = q_initial * ones(N_q, 1);
S_gen_uniform = objective_entropy(q_uniform) * dt;

T_optimal = [T1_P8; T1_P8 + (dt/c) * cumsum(q_optimal(1:end-1))];
T_uniform = [T1_P8; T1_P8 + (dt/c) * cumsum(q_uniform(1:end-1))];

fprintf('P8: S_{gen, optimal}: %.2f kJ/K. S_{gen, uniform}: %.2f kJ/K.\n', S_gen_optimal, S_gen_uniform);

function S_gen_sum = objective_entropy_P8(q, T_A, T_b, c, dt)
    T_curr = T_A;
    S_gen_sum = 0;
    for i = 1:length(q)
        S_gen_dot = (q(i)/T_curr) - (q(i)/T_b);
        S_gen_sum = S_gen_sum + S_gen_dot;
        T_curr = T_curr + (q(i)/c) * dt;
    end
end

figure(5);
t_plot = linspace(0, tf, N_q);
subplot(2,1,1); stairs(t_plot, q_optimal, 'LineWidth', 2); hold on; stairs(t_plot, q_uniform, '--');
title('P8: Optimal vs Uniform Heating Profile'); xlabel('Time (s)'); ylabel('Heating Rate q (kW)'); legend('Optimal q', 'Uniform q'); grid on; hold off;
subplot(2,1,2); plot(linspace(0, tf, N_q+1), [T1_P8; T_optimal], 'LineWidth', 2); hold on; plot(linspace(0, tf, N_q+1), [T1_P8; T_uniform], '--');
title('P8: Temperature Profile'); xlabel('Time (s)'); ylabel('Temperature T (K)'); legend('Optimal T', 'Uniform T'); grid on; hold off;

fprintf('\n');

%% Problem 9: Reaction Heat and Parameter Estimation Using Synthetic Data

A_true = 1e5; E_true = 45; R_P9 = 8.314e-3; c_P9 = 1.8; T_init_P9 = 300;
t_span_P9 = [0, 4000]; noise_level = 0.01;

ode_func_true = @(t, T) (1/c_P9) * (A_true*exp(-E_true./(R_P9*T)) + 2000*exp(-0.001*t));
[t_synth, T_synth_clean] = ode45(ode_func_true, t_span_P9, T_init_P9);
rng(1); T_synth = T_synth_clean .* (1 + noise_level * randn(size(T_synth_clean)));

model_P9 = @(X, t_data) simulate_temperature_P9(X, t_data, c_P9, R_P9, T_init_P9);
[X_estimated, ~] = lsqcurvefit(model_P9, [5e4, 30], t_synth, T_synth, [0, 0], []);

fprintf('P9: A_estimated: %.2e, E_estimated: %.2f kJ/mol.\n', X_estimated(1), X_estimated(2));

function T_sim = simulate_temperature_P9(X, t_data, c, R, T_init)
    ode_func = @(t, T) (1/c) * (X(1)*exp(-X(2)./(R*T)) + 2000*exp(-0.001*t));
    [t_out, T_out] = ode45(ode_func, t_data, T_init);
    T_sim = interp1(t_out, T_out, t_data, 'linear');
end

fprintf('\n');

%% Problem 10: Entropy-Constrained Design of a Heat Pump Stage

Tc_P10 = 270; Th_P10 = 320; alpha = 0.02; k = 50; S_gen_max = 0.05;

COP = @(rp) (Tc_P10 / (Th_P10 - Tc_P10)) * (1 - alpha * (rp - 1).^2 ./ rp);
Wc = @(rp) k * (rp.^0.5 - 1);
S_gen = @(rp) (COP(rp)./Th_P10 - (COP(rp) - 1)./Tc_P10) .* Wc(rp);

objective_P10 = @(rp) -COP(rp);
nonlcon_P10 = @(rp) deal(S_gen(rp) - S_gen_max, []);

rp_min = 1.1; rp_max = 6; rp_guess = 3;
options_P10 = optimoptions('fmincon', 'Display', 'off');
[rp_optimal, max_COP_neg] = fmincon(objective_P10, rp_guess, [], [], [], [], rp_min, rp_max, nonlcon_P10, options_P10);
COP_max = -max_COP_neg;

fprintf('P10: Optimal r_p: %.4f. Max COP: %.4f. S-gen: %.4f kJ/K.\n', rp_optimal, COP_max, S_gen(rp_optimal));

rp_plot = linspace(rp_min, rp_max, 100);
figure(6);
subplot(2,1,1); plot(rp_plot, COP(rp_plot)); hold on; plot(rp_optimal, COP_max, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
title('P10: COP vs Pressure Ratio'); xlabel('Pressure Ratio r_p'); ylabel('COP'); grid on; hold off;
subplot(2,1,2); plot(rp_plot, S_gen(rp_plot)); hold on; yline(S_gen_max, 'r--', 'LineWidth', 2, 'DisplayName', 'Max S-gen');
plot(rp_optimal, S_gen(rp_optimal), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
title('P10: Entropy Generation vs Pressure Ratio'); xlabel('Pressure Ratio r_p'); ylabel('S-gen (kJ/K)'); legend('S-gen(r_p)', 'Max S-gen', 'Location', 'NorthWest'); grid on; hold off;

fprintf('\n');