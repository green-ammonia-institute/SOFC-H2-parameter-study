%% Load current and power data
nSim = 6;
VV_sim = zeros(nSim, 12);
for i = 1:nSim
    tab_current = readtable("Current - Simulation " + i + ".csv");
    tab_power = readtable("Power - Simulation " + i + ".csv");
    VV_sim(i,:) = tab_power{:,1};
    II_sim(i,:) = tab_current{:,2};
    PP_sim(i,:) = tab_power{:,2};
end

%Vs = round(1e2* unique(tab_power{:,1}) ) * 1e-2;

% Load experimental data
tab_exp = readtable("Test_850_H2puro_d.csv");
% Round to 2nd decimal
%V_exp = flip(round(1e2*tab_exp{:,1}) * 1e-2);
V_exp = tab_exp{end:-2:1, 1}.';
% Only use simulated voltages
%[V_exp, ind] = intersect(V_exp, V_sim, "stable");
I_exp = tab_exp{end:-2:1,3}.';
%I_exp = flip(tab_exp{ind,3});
P_exp = tab_exp{end:-2:1, 5}.';
%P_exp = flip(tab_exp{ind,5});

% Set default plot settings
set(0,'defaultAxesFontSize', 18);
set(0,"defaultTextInterpreter", "latex");
set(0,"defaultLegendInterpreter", "latex");

%% Relative change
PP_diff = zeros(nSim-1, 12);
for i = 1:nSim-1
    % Calculate the relative difference curve between simulation i and
    % simulation i+1
    PP_diff(i,:) = abs(PP_sim(i+1,:) - PP_sim(i,:))./PP_sim(i,:);
end
% Maximum relative change
max_diff = max(PP_diff,[],2);

%% Plot runtimes and errors

% Import runtime data
tab_time = readtable("Runtimes.csv");
nElem = tab_time{:,1};
T = tab_time{:,2};
% Colors
left_color = [.5 .5 0];
right_color = [0 .5 .5];
% Default plot settings
set(0,'defaultAxesColorOrder',[left_color; right_color]);
% Plot
fig = figure(Position=[400,200,1000,500]);
hold on;
xlabel("Number of elements");
% Time curve
yyaxis left;
ylabel("Time [s]");
plot(nElem, T, Marker="*");
% Error curve
yyaxis right;
ylabel("Max. relative error");
plot(nElem(2:end), max_diff, Marker="*");
title("Runtime \& relative change vs. number of elements")
% Save plot
ax = gca;
exportgraphics(ax, "Mesh analysis.pdf", BackgroundColor="none", ContentType="vector");

%% Plot polarization and power curves

% figure;
tcl = tiledlayout(1,2);
colors = jet(nSim);
% Polarization
nexttile;
plot(I_exp, V_exp, 'k', 'DisplayName', "Exp");
hold on
for i = 1:nSim
    plot(II_sim(i,:), VV_sim(i,:), Color=colors(i,:), DisplayName="Sim " + i);
end
xlim([0, 2.5]);
ylim([0, 1.3]);
xlabel("Average current [A/cm^2]");
ylabel("Voltage [V]");
legend();
title("Polarization");
% Power
nexttile;
plot(I_exp, P_exp, '-k', 'DisplayName', "Exp");
hold on
for i = 1:nSim
    plot(II_sim(i,:), PP_sim(i,:), Color=colors(i,:), DisplayName="Sim " + i);
end
xlim([0, 2.5]);
ylim([0, 0.7]);
xlabel("Average current [A/cm^2]");
ylabel("Average power [W/cm^2]");
legend();
title("Power");