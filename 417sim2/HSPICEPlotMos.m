% Load results from HSPICE
x = loadsig('MOSSim_HSPICE.sw0');   % Change "netlist" to the appropriate name

% Assign different signals to variables
lssig(x)
% (Use lssig(x) to list possible signals)

% For some reason I can't get the current to show correctly, so I have to
% invert it here (reversing the diode and/or source doesn't work)
%y = -y; 

%Part 2-4

I_D = evalsig(x, 'i_vds');
V_GS = evalsig(x, 'v_2');
I_D = -I_D;

figure(1)
hold on
plot(V_GS,I_D)
plot(0.3007,5.67589e-8,'r*')
legend('V_G_S = 0.3','V_G_S = 0.6','V_G_S = 0.9')
title('Transfer Curves (linear y-axis)');
xlabel('V-GS');
ylabel('I_D');
hold off

%Estimate SS, ss = 1/slope
slope1 = (9.41748e-9-2.0179e-10)/(0.4752-0.4648);
ss1 = 1/slope1;
%2.71 decades between plot so divide by 2.71
ss1 = ss1/1.67;



figure(2)
semilogy(V_GS,I_D)
legend('V_G_S = 0.3, Subthreshold Swing = 6.75e5','V_G_S = 0.6','V_G_S = 0.9')
title('Transfer Curves (log scale y-axis)');
xlabel('V-GS');
ylabel('I_D');


