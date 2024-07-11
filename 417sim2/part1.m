%ECE 417 Sim 2 - MATLAB
%Constants
G_ox = 2.5e-6;
Q_f = 1 * 10^11; %cm^-2
T = 300; %Kelvin
q = 1.60217663e-19; %Coulombs
K_B = 1.3806452e-23;
N_A = 9e17; %cm-3 for substrate
eps_0 = 8.854187817e-12;
eps_0x = 3.9;
n_i = 1e10; %cm^-3 for substrate
V_BS = 0;
E_g = 1.12; %eV
eps_si = 11.9;
V_GS = 0;
L = 0.5e-6; % channel length;
W = 0.8e-6; % channel width;
%Part 1
%Estimate V_TO
C_ox = (eps_0 * eps_0x)/ G_ox;
gamma = (sqrt(2*q*eps_si*eps_0*N_A)/C_ox);
Psi_B = ((K_B*T)/q)*log(N_A/n_i);
Phi_MS = -1*(E_g/2) + -1*Psi_B;
V_FB = Phi_MS - ((q*Q_f)/C_ox);
V_TO = V_FB + (2*Psi_B) + (gamma*(sqrt((2*Psi_B)-V_BS)));
%Part 2m
%Output Curves
%Calculate V_dsat
%Change for given V_GS
%K = (sqrt(q*eps_si*N_A))/C_ox;
K = gamma/(sqrt(2));
M = 1 + (K/(sqrt(2*(2*Psi_B-V_BS))));
V_dsat = (V_GS - V_TO) / M;
%As V_GS changes, V_dsat will change for each iteration of V_GS
%The sweep of V_DS changes the region of operation for a given V_dsat
% sweep V_DS from 0V to 2V
    %V_GS > V_T0 and V_DS < V_dsat
   %I_DS_pre_sat = (W/L)*mu*C_ox((V_GS-V_TO)*V_DS - (M*(V_DS^2)/L));
   %V_GS >= V_T0 and V_DS >= V_dsat
   %I_DS_sat = (W/(2*M*L))*mu*C_ox*((V_GS-V_TO)^2);
% Constants and parameters
W_over_L = W/L; % Width-to-Length ratio
mu = 0.1; % Mobility
% Voltage range
V_DS_range = 0:0.01:2; % VDS range from 0 to 2 V
% Choose a specific value for V_GS
V_GS1 = 0;
% Calculate V_DSat
V_DSat1 = (V_GS1 - V_TO) / M;
% Cutoff region
I_D_cutoff1 = zeros(size(V_DS_range));
I_D_cutoff1(V_GS1 < V_TO) = 0;
% Pre-saturation region
I_DS_pre_sat1 = (W_over_L) * mu * C_ox * ((V_GS1 - V_TO) * V_DS_range - (M * (V_DS_range.^2) / 2));
%I_DS_pre_sat1(V_DS_range < V_DSat1) = (W_over_L) * mu * C_ox * ((V_GS1 - V_TO) * V_DS_range - (M * (V_DS_range.^2) / 2));
% Saturation region
I_DS_sat1 = (W_over_L / (2 * M)) * mu * C_ox * ((V_GS1 - V_TO)^2);
I_DS_sat1(V_DS_range >= V_DSat1) = (W_over_L / (2 * M)) * mu * C_ox * ((V_GS1 - V_TO)^2);
% Choose a specific value for V_GS
V_GS2 = 1.3;
% Calculate V_DSat
V_DSat2 = (V_GS2 - V_TO) / M;
% Cutoff region
I_D_cutoff2 = zeros(size(V_DS_range));
I_D_cutoff2(V_GS2 < V_TO) = 0;
% Pre-saturation region
I_DS_pre_sat2 = (W_over_L) * mu * C_ox * ((V_GS2 - V_TO) * V_DS_range - (M * (V_DS_range.^2) / 2));
%I_DS_pre_sat2(V_DS_range < V_DSat2) = (W_over_L) * mu * C_ox * ((V_GS2 - V_TO) * V_DS_range - (M * (V_DS_range.^2) / 2));
% Saturation region
I_DS_sat2 = (W_over_L / (2 * M)) * mu * C_ox * ((V_GS2 - V_TO)^2);
I_DS_sat2(V_DS_range >= V_DSat2) = (W_over_L / (2 * M)) * mu * C_ox * ((V_GS2 - V_TO)^2);
% Choose a specific value for V_GS
V_GS3 = 1.6;
% Calculate V_DSat
V_DSat3 = (V_GS3 - V_TO) / M;
% Cutoff region
I_D_cutoff3 = zeros(size(V_DS_range));
I_D_cutoff3(V_GS3 < V_TO) = 0;
% Pre-saturation region
I_DS_pre_sat3 = (W_over_L) * mu * C_ox * ((V_GS3 - V_TO) * V_DS_range - (M * (V_DS_range.^2) / 2));
%I_DS_pre_sat3(V_DS_range < V_DSat3) = (W_over_L) * mu * C_ox * ((V_GS3 - V_TO) * V_DS_range - (M * (V_DS_range.^2) / 2));
% Saturation region
I_DS_sat3 = (W_over_L / (2 * M)) * mu * C_ox * ((V_GS3 - V_TO)^2);
I_DS_sat3(V_DS_range >= V_DSat3) = (W_over_L / (2 * M)) * mu * C_ox * ((V_GS3 - V_TO)^2);
% Choose a specific value for V_GS
V_GS4 = 2.2;
% Calculate V_DSat
V_DSat4 = (V_GS4 - V_TO) / M;
% Cutoff region
I_D_cutoff4 = zeros(size(V_DS_range));
I_D_cutoff4(V_GS4 < V_TO) = 0;
% Pre-saturation region
I_DS_pre_sat4 = (W_over_L) * mu * C_ox * ((V_GS4 - V_TO) * V_DS_range - (M * (V_DS_range.^2) / 2));
%I_DS_pre_sat4(V_DS_range < V_DSat4) = (W_over_L) * mu * C_ox * ((V_GS4 - V_TO) * V_DS_range - (M * (V_DS_range.^2) / 2));
% Saturation region
I_DS_sat4 = (W_over_L / (2 * M)) * mu * C_ox * ((V_GS4 - V_TO)^2);
I_DS_sat4(V_DS_range >= V_DSat4) = (W_over_L / (2 * M)) * mu * C_ox * ((V_GS4 - V_TO)^2);
% Choose a specific value for V_GS
V_GS5 = 2.9;
% Calculate V_DSat
V_DSat5 = (V_GS5 - V_TO) / M;
% Cutoff region
I_D_cutoff5 = zeros(size(V_DS_range));
I_D_cutoff5(V_GS5 < V_TO) = 0;
% Pre-saturation region
I_DS_pre_sat5 = (W_over_L) * mu * C_ox * ((V_GS5 - V_TO) * V_DS_range - (M * (V_DS_range.^2) / 2));
%I_DS_pre_sat5(V_DS_range < V_DSat5) = (W_over_L) * mu * C_ox * ((V_GS5 - V_TO) * V_DS_range - (M * (V_DS_range.^2) / 2));
% Saturation region
I_DS_sat5 = (W_over_L / (2 * M)) * mu * C_ox * ((V_GS5 - V_TO)^2);
I_DS_sat5(V_DS_range >= V_DSat5) = (W_over_L / (2 * M)) * mu * C_ox * ((V_GS5 - V_TO)^2);
% Plot the curves
figure;
%plot(V_DS_range(V_DS_range < V_DSat1), I_D_cutoff1(V_DS_range < V_DSat1));
%plot(V_DS_range(V_DS_range < V_DSat2), I_D_cutoff2(V_DS_range < V_DSat2));
%plot(V_DS_range(V_DS_range < V_DSat3), I_D_cutoff3(V_DS_range < V_DSat3));
%plot(V_DS_range(V_DS_range < V_DSat4), I_D_cutoff4(V_DS_range < V_DSat4));
%plot(V_DS_range(V_DS_range < V_DSat5), I_D_cutoff5(V_DS_range < V_DSat5));
hold on;
plot(V_DS_range(V_DS_range < V_DSat1), I_DS_pre_sat1(V_DS_range < V_DSat1), 'DisplayName', '0V');
plot(V_DS_range(V_DS_range < V_DSat2), I_DS_pre_sat2(V_DS_range < V_DSat2), 'DisplayName', '1.3V pre sat');
plot(V_DS_range(V_DS_range < V_DSat3), I_DS_pre_sat3(V_DS_range < V_DSat3), 'DisplayName', '1.6V pre sat');
plot(V_DS_range(V_DS_range < V_DSat4), I_DS_pre_sat4(V_DS_range < V_DSat4), 'DisplayName', '2.2V pre sat');
plot(V_DS_range(V_DS_range < V_DSat5), I_DS_pre_sat5(V_DS_range < V_DSat5), 'DisplayName', '2.9V pre sat');
plot(V_DS_range(V_DS_range >= V_DSat1), I_DS_sat1(V_DS_range >= V_DSat1), 'DisplayName', '0V');
plot(V_DS_range(V_DS_range >= V_DSat2), I_DS_sat2(V_DS_range >= V_DSat2), 'DisplayName', '1.3V sat');
plot(V_DS_range(V_DS_range >= V_DSat3), I_DS_sat3(V_DS_range >= V_DSat3), 'DisplayName', '1.6V sat');
plot(V_DS_range(V_DS_range >= V_DSat4), I_DS_sat4(V_DS_range >= V_DSat4), 'DisplayName', '2.2V sat');
plot(V_DS_range(V_DS_range >= V_DSat5), I_DS_sat5(V_DS_range >= V_DSat5), 'DisplayName', '2.9V sat');
% Customize the plot
title(sprintf('NMOS ID – VDS Curve for VGS = 0,1.3,1.6,2.2,2.9 V'));
xlabel('VDS (Volts)');
ylabel('ID (Amperes)');
legend('Location', 'Best');
grid on;
hold off;
%Plot transfer curves
%Same as before except static V_DS and sweep of V_GS
%For V_GS < V_T0
%I_D = 0;
%For V_GS >= V_T0 and V_DS < V_dsat
%I_DS = (W/L)*mu*C_ox;
%For V_GS >= V_T0 and V_DS >= V_dsat
%I_DS = (W/(2*M*L))*mu*C_ox*((V_GS-V_TO)^2);
%Plot
% Voltage ranges
V_GS_range = 0:0.01:3; % VGS range from 0 to 3 V
V_DS_values = [0.3, 0.6, 0.9]; % VDS values
V_BS = 0; % Back-gate voltage
% Initialize arrays to store data for sub-threshold swing calculation
VGS_subthreshold = [];
IDS_subthreshold = [];
figure;
% Linear y-axis scale plot
subplot(1, 2, 1);
for i = 1:length(V_DS_values)
   V_DS = V_DS_values(i);
  
   % Calculate ID for the given VGS and VDS
   ID = (W_over_L) * mu * C_ox * ((V_GS_range - V_TO) .* V_DS - (M * (V_DS^2) / 2));
  
   % Plot the transfer curve
   plot(V_GS_range, ID, 'DisplayName', sprintf('VDS = %.1f V', V_DS));
   hold on;
  
   % Store data for sub-threshold swing calculation (for VDS = 0.3 V)
   if V_DS == 0.3
       VGS_subthreshold = [VGS_subthreshold, V_GS_range];
       IDS_subthreshold = [IDS_subthreshold, ID];
   end
end
title('ID – VGS Transfer Curve (Linear Scale)');
xlabel('VGS (Volts)');
ylabel('ID (Amperes)');
legend('Location', 'Best');
grid on;
% Log y-axis scale plot
subplot(1, 2, 2);
% Choose VDS = 0.3 V for sub-threshold swing calculation
V_DS_subthreshold = 0.3;
ID_subthreshold = (W_over_L) * mu * C_ox * ((V_GS_range - V_TO) .* V_DS_subthreshold - (M * (V_DS_subthreshold^2) / 2));
V_DS_subthreshold1 = 0.6;
ID_subthreshold1 = (W_over_L) * mu * C_ox * ((V_GS_range - V_TO) .* V_DS_subthreshold1 - (M * (V_DS_subthreshold1^2) / 2));
V_DS_subthreshold2 = 0.9;
ID_subthreshold2 = (W_over_L) * mu * C_ox * ((V_GS_range - V_TO) .* V_DS_subthreshold2 - (M * (V_DS_subthreshold2^2) / 2));
semilogy(V_GS_range, ID_subthreshold, 'DisplayName', sprintf('VDS = %.1f V', V_DS_subthreshold));
semilogy(V_GS_range, ID_subthreshold1, 'DisplayName', sprintf('VDS = %.1f V', V_DS_subthreshold1));
semilogy(V_GS_range, ID_subthreshold2, 'DisplayName', sprintf('VDS = %.1f V', V_DS_subthreshold2));
title('ID – VGS Transfer Curve (Log Scale)');
xlabel('VGS (Volts)');
ylabel('ID (Amperes)');
legend('Location', 'Best');
grid on;
hold off;
% (a) Estimate VT from the linear y-axis scale transfer curve
VT_estimate = V_GS_range(find(ID_subthreshold > 1e-12, 1, 'first'));
% (b) Estimate sub-threshold swing (in mV/dec) from the log y-axis scale
subthreshold_indices = ID_subthreshold > 1e-12;
subthreshold_slope = diff(log(ID_subthreshold(subthreshold_indices))) / diff(V_GS_range(subthreshold_indices));
subthreshold_swing_mV_per_dec = -1e3 / subthreshold_slope;
fprintf('(a) VT (Threshold Voltage) Estimate: %.3f V\n', VT_estimate);
fprintf('(b) Sub-threshold Swing Estimate: %.2f mV/dec\n', subthreshold_swing_mV_per_dec);
