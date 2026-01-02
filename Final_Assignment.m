% Given Data
T1 = 298;
p1 = 1e5;
eta_c = 0.88;
eta_gt = 0.90;
eta_st = 0.92;
eta_p = 0.60;

T_cond = 303;
dT_pinch = 10;
x_min = 0.92;

LHV = 50e6;
cp = 1005;
gamma = 1.4;

T3_range = 900:100:1400;

eta_comb = zeros(size(T3_range));
eta_br = zeros(size(T3_range));
eta_ra = zeros(size(T3_range));
eta_rp = zeros(size(T3_range));
eta_p7 = zeros(size(T3_range));

%% Main Optimization Loop
for i = 1:length(T3_range)
    T3 = T3_range(i) + 273.15;
    best_eta = -inf;

    for rp = 5:0.5:30
        for p7 = 5e5:5e5:4e6

            [eta_cyc, eta_b, eta_r, valid, states] = combined_cycle(T1,p1,T3,rp,p7,eta_c,eta_gt,eta_st,eta_p,T_cond,dT_pinch,x_min);

            if valid && eta_cyc > best_eta
                best_eta = eta_cyc;
                eta_br(i) = eta_b;
                eta_ra(i) = eta_r;
                opt_rp(i) = rp;
                opt_p7(i) = p7;
                best_states = states;
            end
        end
    end

    eta_comb(i) = best_eta;
end

Results = table( ...
    T3_range', ...
    opt_rp', ...
    eta_br'*100, ...
    eta_ra'*100, ...
    eta_comb'*100, ...
    'VariableNames', {'T3_C','Optimal_rp','Eta_Brayton(%)','Eta_Rankine(%)','Eta_Combined(%)'});

disp('OPTIMIZATION RESULTS')
disp(Results)

%% Ploting the results
figure;
plot(T3_range, eta_comb*100, 'LineWidth',2)
xlabel('T_3 (°C)')
ylabel('Combined Cycle Efficiency (%)')
title('Combined Cycle Efficiency vs Turbine Inlet Temperature')
grid on

% T-s Diagram
figure;
plot(best_states.s_brayton, best_states.T_brayton, '-o', 'LineWidth',2)
hold on
plot(best_states.s_rankine, best_states.T_rankine, '-o', 'LineWidth',2)

for i = 1:length(best_states.s_brayton)
    text(best_states.s_brayton(i), best_states.T_brayton(i), [' ', num2str(i)], 'FontSize',10);
end

for i = 1:length(best_states.s_rankine)
    text(best_states.s_rankine(i), best_states.T_rankine(i), [' ', num2str(i+4)], 'FontSize',10);
end

xlabel('Entropy (kJ/kgK)')
ylabel('Temperature (K)')
title('T–s Diagram (Combined Cycle)')
legend('Brayton Cycle','Rankine Cycle','Location','best')
grid on

% p-v Diagram
figure;
plot(best_states.v_brayton, best_states.p_brayton/1e5,'-o','LineWidth',2)
hold on
plot(best_states.v_rankine, best_states.p_rankine/1e5,'-o','LineWidth',2)

for i = 1:length(best_states.v_brayton)
    text(best_states.v_brayton(i), best_states.p_brayton(i)/1e5, ['  ', num2str(i)], 'FontSize',10);
end

for i = 1:length(best_states.v_rankine)
    text(best_states.v_rankine(i), best_states.p_rankine(i)/1e5, ['  ', num2str(i+4)], 'FontSize',10);
end

xlabel('Specific Volume (m^3/kg)')
ylabel('Pressure (bar)')
title('p–v Diagram (Combined Cycle)')
legend('Brayton Cycle','Rankine Cycle','Location','best')
grid on
           