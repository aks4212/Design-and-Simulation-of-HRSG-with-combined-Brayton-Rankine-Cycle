function [eta_comb, eta_br, eta_ra, valid, states] = combined_cycle(T1,p1,T3,rp,p7,eta_c,eta_gt,eta_st,eta_p,T_cond,dT_pinch,x_min)

valid = true;

cp = 1005;
gamma = 1.4;
R = 287;

%% Compressor
p2 = p1 * rp;
T2s = T1*(p2/p1)^((gamma-1)/gamma);
T2 = T1 + (T2s - T1)/eta_c;

Wc = cp*(T2 - T1);

%% Combustor
Qin = cp*(T3 - T2);

%% Gas Turbine
T4s = T3*(1/rp)^((gamma-1)/gamma);
T4 = T3 - eta_gt*(T3 - T4s);
Wt = cp*(T3 - T4);

eta_br = (Wt - Wc)/Qin;

%% HRSG
T_gas_exit = T4 - dT_pinch;

% steam-tables data
hf = 419e3; 
hfg = 2257e3;

T10 = T_gas_exit - dT_pinch;
h10 = hf + hfg;

Q_available = cp*(T4 - T_gas_exit);
Q_required = h10 - hf;
steam_ratio = Q_available / Q_required;

x = x_min;
h11 = hf + x*hfg;

%% Steam Turbine
Wt_s = steam_ratio*eta_st*(h10 - h11);

%% Pump
Wp = (p7 - p1)/1000/eta_p;

%% Rankine Efficiency
eta_ra = (Wt_s - Wp)/(steam_ratio*(h10 - hf));

%% Total Efficiency
eta_comb = (Wt + Wt_s - Wc - Wp)/Qin;


if x < x_min || eta_comb < 0
    valid = false;
end

% for plots--
states.T_brayton = [T1 T2 T3 T4];
states.s_brayton = [1 2 3 4];     
states.p_brayton = [p1 p2 p2 p1];
states.v_brayton = R*states.T_brayton./states.p_brayton;

states.T_rankine = [T_cond T10 T_cond];
states.s_rankine = [5 6 7];
states.p_rankine = [p1 p7 p1];
states.v_rankine = R*states.T_rankine./states.p_rankine;

end