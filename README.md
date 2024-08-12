# Moisture-on-Porous-Media-Project
---

### Step 1 **Setup and Constants Definition**
   - **Gas Constants**: Calculate specific gas constants for vapor (`R_v`) and air (`R_a`) using the universal gas constant `B` and respective molar masses (`M_v` for vapor and `M_a` for air).
   - **Density Calculation (`Rho_m`)**: Define the mixture gas constant `R_m` based on the vapor mass fraction `Y_v_val`, and calculate the overall density of the mixture `Rho_m` using the ideal gas law.

   ```matlab
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
   ```

### Step 2 **Boundary Conditions for Temperature and Humidity**
   - **Initial and Boundary Conditions**: Define the temperature (`T_in`, `T_out`) and relative humidity (`RH_in`, `RH_out`) at the inner and outer surfaces of the wall, then calculate the vapor pressure, vapor density, and air density at both boundaries.
   - **Vapor Mass Fraction (`Y_in` and `Y_out`)**: Calculate the vapor mass fraction at the inner and outer surfaces.

   ```matlab
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
   ```

### Step 3 **Temperature and Diffusion Model**
   - **Temperature Profile (`T_wall`)**: Define the linear temperature profile across the wall thickness.
   - **Diffusion Equation**: Use Fick's law to define the vapor flux `j_v` and set up the differential equation for the vapor mass fraction `Y_v_val` within the wall.
   - **Solution of the Diffusion Equation**: Solve the diffusion equation under specified boundary conditions to get the vapor mass fraction profile `Y_v_val(x)`.

   ```matlab
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
   ```

### Step 4 **Relative Humidity Profile**
   - **Saturation Vapor Pressure**: Calculate the saturation vapor pressure at the wall temperature.
   - **Relative Humidity (`RH_wall`)**: Define the relative humidity inside the wall as a function of position.

   ```matlab
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
   ```

### Step 5 **Moisture Content using GAB Equation**
   - **GAB Equation**: Compute the moisture content within the wall using the GAB (Guggenheim-Anderson-de Boer) equation.

   ```matlab
   wm = 0.028;
   C = 10;
   K = 0.62;
   w_GAB = simplify((wm*C*K.*RH_wall) / ((1-K.*RH_wall) * (1+K*(C-1)*RH_wall)),...
       "Seconds",240);

   figure
   fplot(w_GAB, [0 L])
   title('Moisture Content (kg/kg)')
   xlabel("wall thickness (m)")
   ```

### Step 6 **Heat Transfer Analysis**
   - **Thermal Conductivity Variation**: Model the thermal conductivity of the wall material as a function of moisture content.
   - **Heat Flux (`q`)**: Define the heat flux within the wall and solve the steady-state heat conduction equation to obtain the temperature distribution.
   - **Temperature Adjustment**: Adjust the temperature profile based on the heat flux.

   ```matlab
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
   ```

### Step 7 **3D Visualization**
   - **Variation of Vapor Mass Fraction with Temperature**: Simulate and visualize the variation of vapor mass fraction with the external temperature across the wall.
   - **Variation of Vapor Mass Fraction with Relative Humidity**: Simulate and visualize the vapor mass fraction with varying external relative humidity.

   ```matlab
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
       Rho_a_out = P_atm/(

R_a*T_out);
       Y_out = Rho_v_out/(Rho_a_out+Rho_v_out);

       T_wall = ((-T_in + T_out) / L) * x + T_in;
       R_m = Y_v_val * R_v + (1 - Y_v_val) * R_a;
       Rho_m = P_atm / (R_m * T_wall);

       j_v = -Rho_m * D_wall * diff(Y_v_val, x);
       Y_v_val_3D = simplify(dsolve(EQ_Y, BC_Y));

       Y_v_matrix(:, i) = double(subs(Y_v_val_3D, x, x_vals));
   end

   figure
   surf(x_vals, T_out_vals, Y_v_matrix', 'EdgeColor', 'none');
   title('Variation of Y_v with T_out across wall thickness')
   xlabel("wall thickness (m)")
   ylabel("T out (K)")
   zlabel("Y_v")

   RH_out_vals = linspace(0.1, 0.9, 20);
   Y_v_matrix_RH = zeros(50, length(RH_out_vals));

   for i = 1:length(RH_out_vals)
       RH_out = RH_out_vals(i);

       P_sat_out3 = P_sat_3D(T_out);
       P_v_out = P_sat_out3 * RH_out;
       Rho_v_out = P_v_out/(R_v*T_out);
       Rho_a_out = P_atm/(R_a*T_out);
       Y_out = Rho_v_out/(Rho_a_out+Rho_v_out);

       Y_v_val_3D = simplify(dsolve(EQ_Y, BC_Y));

       Y_v_matrix_RH(:, i) = double(subs(Y_v_val_3D, x, x_vals));
   end

   figure
   surf(x_vals, RH_out_vals, Y_v_matrix_RH', 'EdgeColor', 'none');
   title('Variation of Y_v with RH_out across wall thickness')
   xlabel("wall thickness (m)")
   ylabel("RH out")
   zlabel("Y_v")
   ```

### Step 8 **Plotting**
   - **Multiple Figures**: The code includes multiple plotting sections that visualize the profiles of vapor mass fraction, relative humidity, moisture content, thermal conductivity, and temperature distribution within the wall.

   ```matlab
   % No additional code needed here, as the plotting is already integrated within each section above.
   ```

---
