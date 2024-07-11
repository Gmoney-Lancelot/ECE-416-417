%Part 2-3 - Output Curves


y = evalsig(x, 'v_1');
w = evalsig(x, 'i_vds');
w = -w;
V = evalsig(x, 'VOLTS');

%Plot results
%V_T0 = 0.490 - from part 1
%M = 1.2051 - from part 1
%Mark plots for cut-off, pre-sat, and sat
%V_DSAT =(V_GS - V_T0)/M
V_T0 = 0.3007;
M = 1.2051;
%V_DSAT0 = (0 - V_T0)/M;
V_DSAT13 =(1.3 - V_T0)/M;
V_DSAT16 =(1.6 - V_T0)/M;
V_DSAT19 =(1.9 - V_T0)/M;
V_DSAT22 =(2.2 - V_T0)/M;

figure
hold on
plot(V,w)
plot(0.67,4.50263e-5,'r*')
plot(0.92,8.3697e-5, 'r*')
plot(1.17,0.000134158, 'r*')
plot(1.42,0.000196777, 'r*')

legend('V_g_s = 0','V_g_s = 1.3','V_g_s = 1.6','V_g_s = 1.9','V_g_s = 2.2')
title('Output Curves');
xlabel('V-DS');
ylabel('I_D');
hold off