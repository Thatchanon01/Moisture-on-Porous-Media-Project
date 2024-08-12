clear;

syms Y_v(x) Y_v_val(x) T x hm

B = 8314;
M_v = 18;
M_a = 28.9;

R_v = B/M_v;
R_a = B/M_a;

R_m = Y_v_val*R_v +(1-Y_v_val)*R_a;
P_atm = 101325; 
Rho_m = P_atm/(R_m*T);
L = 0.1; % Assume wall thickness

P_sat = 611 * exp(17.08*(T-273.150)/(T-38.97));

T_out = 25+273.15;
RH_out = 0.5;
P_sat_out = subs(P_sat,T_out); %[kPa]
P_v_out = P_sat_out * RH_out;
Rho_v_out = P_v_out/(R_v*T_out);
Rho_a_out = P_atm/(R_a*T_out);
Y_out = Rho_v_out/(Rho_a_out+Rho_v_out);
Y_out_val = double(Y_out);

T_in = 22+273.15; %Assumption
RH_in = 0.6; %Assumption
P_sat_in = subs(P_sat,T_in);
P_v_in = P_sat_in * RH_in;
Rho_v_in = P_v_in/(R_v*T_in);
Rho_a_in = (P_atm-P_v_in)/(R_a*T_in);
Y_in = Rho_v_in/(Rho_a_in+Rho_v_in);
Y_in_val = double(Y_in);

Rho_a = 1.208;
cp_a = 1.007;
cp_w = 1.86;
h_v = 2501;
Alpha = 2.18e-5; %thermal diffusity
D_av = 1.87e-10*((T^2.072)*(P_atm/101325)); %Air-Vapour P[atm] T[K]

T_wall = ((-T_in + T_out) / L) * x + T_in;
D_wall = 2.425e-12;  % D vapour-mat
Rho_m = subs(Rho_m,T,T_wall);

j_v = -Rho_m*D_wall*diff(Y_v_val, x);

BC_Y = [Y_v_val(0) == Y_in_val; Y_v_val(L) == Y_out_val];
EQ_Y = 0 == diff(-j_v,x);
Y_v_val = simplify(dsolve(EQ_Y, BC_Y),"Seconds",180);

figure
fplot(Y_v_val, [0 L])
title('Vapour mass fraction')
xlabel("wall thickness (m)")

Y_wall_in = double(subs(Y_v_val, x, 0));
Y_wall_out = double(subs(Y_v_val, x, L));

%% RH wall & GAB equation

P_sat_wall = subs(P_sat,T_wall); %[kPa]
 
Y_v_subs = Y_v_val;
R_m_subs = Y_v_subs*R_v +(1-Y_v_subs)*R_a;
Rho_m_subs = P_atm /(R_m_subs*T_wall);

Rho_v = Rho_m_subs*Y_v_subs;

P_v = simplify(Rho_v*R_v*T_wall, "Seconds",120); 

RH_wall = simplify(P_v / P_sat_wall,"Seconds",180);

figure
fplot(RH_wall, [0 L])
title('RH_wall profile [%]')
xlabel("wall thickness (m)")

wm = 0.028;
C = 10;
K = 0.62;
w_GAB = simplify((wm*C*K.*RH_wall) / ((1-K.*RH_wall) * (1+K*(C-1)*RH_wall)),...
    "Seconds",240);

figure
fplot(w_GAB, [0 L])
title('Moisture Content (kg/kg)')
xlabel("wall thickness (m)")

%% Heat transfer
%Case k varies based on Moisture content and density

syms T(x)

b = 4; %thermal conductivity supplement
Rho_s = 416; %Assume material
k_0 = 0.125; %thermal conductivity
w_vol = w_GAB*Rho_s; %[kg/m3]
k_w = k_0 * (1 + b * w_vol / Rho_s);

k_perc = (double(subs(k_w,x,L/2)) - k_0)*100;
q = -k_w * diff(T,x);
EQ_T = 0 == diff(-q,x); %steady state

x_val = linspace(0, L, 11);

q_wall = -k_w * diff(T_wall,x);
q_subs = double(subs(q_wall,x,x_val));

figure
fplot(k_w, [0 L])
title('Thermal Conductivity (W/mK)')
xlabel("wall thickness (m)")

HeatTransfer = q_subs(:, 6);

BC_T = [T(0) == T_in T(L) == T_out];

dT = q_subs.*x_val./k_w;
T_after = subs(T_wall - dT, x, 0);

figure
plot(x_val, T_after)
title('Temperature_wall Profile (K)')
xlabel("wall thickness (m)")


%% 3D Visualization

P_sat_3D = @(T) 611 * exp(17.08 * (T - 273.15) / (T - 38.97));

P_sat_in_3D = P_sat_3D(T_in);
P_v_in = P_sat_in_3D * RH_in;
Y_in = P_v_in / (P_atm - P_v_in);

T_out_vals = linspace(-5 + 273.15, 37 + 273.15, 20); 
Y_v_matrix = zeros(50, length(T_out_vals)); 

x_vals = linspace(0, L, 50);

for i = 1:length(T_out_vals)
    T_out = T_out_vals(i);

    P_sat_out3 = P_sat_3D(T_out);
    P_v_out = P_sat_out3 * 0.5; % Assume RH_out = 0.5
    Rho_v_out = P_v_out/(R_v*T_out);
    Rho_a_out = P_atm/(R_a*T_out);
    Y_out = Rho_v_out/(Rho_a_out+Rho_v_out);

    T_wall = ((-T_in + T_out) / L) * x + T_in;

    R = Y_v * R_v + (1 - Y_v) * R_a;
    Rho = P_atm / (R * T_wall);
    D_av = 1.87e-10 * ((T_wall^2.072) * (P_atm / 101325));
    j_v = -Rho * D_av * diff(Y_v, x);

    BC_Y = [Y_v(0) == Y_in; Y_v(L) == Y_out];
    EQ_Y = diff(j_v, x) == 0;
    Y_v_sol = dsolve(EQ_Y, BC_Y);
    Y_v_sol = simplify(Y_v_sol);

    for j = 1:length(x_vals)
        Y_v_matrix(j, i) = double(subs(Y_v_sol, x, x_vals(j)));
    end
end

[X, T_OUT] = meshgrid(x_vals, T_out_vals);
Y_V = Y_v_matrix';

figure;
surf(X, T_OUT, Y_V);
xlabel('Position (m)');
ylabel('Outside Temperature (K)');
zlabel('Vapour Mass Fraction (Y_v)');
title('Vapour Mass Fraction (Y_v) Variation with Position and Outside Temperature');
grid on;

x_positions = [0.025, 0.05, 0.075];
colors = ['r', 'g', 'b'];

figure;
hold on;
for k = 1:length(x_positions)
    [~, idx] = min(abs(x_vals - x_positions(k)));
    plot(T_out_vals, Y_V(:, idx), colors(k));
end
xlabel('Outside Temperature (K)');
ylabel('Vapour Mass Fraction (Y_v)');
legend('x = 0.025', 'x = 0.05', 'x = 0.075');
title('Vapour Mass Fraction at Different Positions with Varying T');
hold off;

RH_out_vals = linspace(0.3, 1, 20); 

Y_v_matrix_RH = zeros(50, length(RH_out_vals));

for i = 1:length(RH_out_vals)
    RH_out = RH_out_vals(i);
    T_out = 25 + 273.15; 

    P_sat_out3 = P_sat_3D(T_out);
    P_v_out = P_sat_out3 * RH_out;
    Rho_v_out = P_v_out/(R_v*T_out);
    Rho_a_out = P_atm/(R_a*T_out);
    Y_out = Rho_v_out/(Rho_a_out+Rho_v_out);

    T_wall = ((-T_in + T_out) / L) * x + T_in;

    R = Y_v * R_v + (1 - Y_v) * R_a;
    Rho = P_atm / (R * T_wall);
    D_av = 1.87e-10 * ((T_wall^2.072) * (P_atm / 101325));
    j_v = -Rho * D_av * diff(Y_v, x);

    BC_Y = [Y_v(0) == Y_in; Y_v(L) == Y_out];
    EQ_Y = diff(j_v, x) == 0;
    Y_v_sol = dsolve(EQ_Y, BC_Y);
    Y_v_sol = simplify(Y_v_sol);

    for j = 1:length(x_vals)
        Y_v_matrix_RH(j, i) = double(subs(Y_v_sol, x, x_vals(j)));
    end
end

[X, RH_OUT] = meshgrid(x_vals, RH_out_vals);
Y_V_RH = Y_v_matrix_RH';

figure;
surf(X, RH_OUT, Y_V_RH);
xlabel('Position (m)');
ylabel('Outside Relative Humidity');
zlabel('Vapour Mass Fraction (Y_v)');
title('Vapour Mass Fraction Variation with Position and Outside Relative Humidity');
grid on;

figure;
hold on;
for k = 1:length(x_positions)
    [~, idx] = min(abs(x_vals - x_positions(k)));
    plot(RH_out_vals, Y_V_RH(:, idx), colors(k));
end
xlabel('Outside Relative Humidity');
ylabel('Vapour Mass Fraction (Y_v)');
legend('x = 0.025', 'x = 0.05', 'x = 0.075');
title('Vapour Mass Fraction at Different Positions with Varying RH');
hold off;
