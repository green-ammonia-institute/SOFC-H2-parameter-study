%% Load current and power data
tab_current = readtable("Current.csv");
tab_power = readtable("Power.csv");
Vs = round(1e2* unique(tab_power{:,1}) ) * 1e-2;

% Set default plot settings
set(0,'defaultAxesFontSize', 18);
set(0,"defaultTextInterpreter", "latex");
set(0,"defaultLegendInterpreter", "latex");

%% Pre-process data
[Xa, Xc, Y_I] = preprocess_data(tab_current);
[Xa, Xc, Y_P] = preprocess_data(tab_power);         % The grids are identical
[nV, nc, na] = size(Y_P);

%% Plot surfaces
for i = 1:nV
    figure;
    tcl = tiledlayout(1,2);
    % Plot current
    nexttile;
    surf(Xa, Xc, squeeze(Y_I(i,:,:)));
    zlim([0, 2.5]);
    xlabel("I_{a,0} [A/cm\textsuperscript{2}]");
    ylabel("I_{c,0} [A/cm\textsuperscript{2}]");
    title("Average cell current density [A/cm\textsuperscript{2}]");
    % Plot power
    nexttile;
    surf(Xa, Xc, squeeze(Y_P(i,:,:)));
    zlim([0, 0.8]);
    xlabel("I_{a,0} [A/cm\textsuperscript{2}]");
    ylabel("I_{c,0} [A/cm\textsuperscript{2}]");
    title("Cell power [W/cm\textsuperscript{2}]");
    % Figure title
    title(tcl, Vs(i) + "V");
end

%% Find exchange currents
% Load experimental data
tab_exp = readtable("Test_850_H2puro_d.csv");
% Round to 2nd decimal
V_exp = flip(round(1e2*tab_exp{:,1}) * 1e-2);
% Only use simulated voltages
[V_exp, ind] = intersect(V_exp, Vs, "stable");
I_exp = flip(tab_exp{ind,3});
P_exp = flip(tab_exp{ind,5});
% Interpolate
Iq = linspace(min(I_exp), max(I_exp), 30);
V_obj = spline(I_exp, V_exp, Iq);
P_obj = spline(I_exp, P_exp, Iq);

% Calculate power errors
err_l1 = zeros(nc, na);
err_l2 = zeros(nc, na);
err_sup = zeros(nc, na);
err_rel = zeros(nc, na);
for j = 1:nc
    for k = 1:na
        P = spline(Y_I(:,j,k), Y_P(:,j,k), Iq);
        err_l1(j,k) = norm(P - P_obj, 1);
        err_l2(j,k) = norm(P - P_obj, 2);
        err_sup(j,k) = max(abs(P - P_obj));
        err_rel(j,k) = mean(abs(P-P_obj)./P_obj);
    end
end

% Find minima
% Objective
min_l1 = min(err_l1, [], "all");
min_l2 = min(err_l2, [], "all");
min_sup = min(err_sup, [], "all");
min_rel = min(err_rel,[], "all");
% Indices must be swapped
[x_l1, y_l1] = find(err_l1 == min_l1);
[x_l2, y_l2] = find(err_l2 == min_l2);
[x_sup, y_sup] = find(err_sup == min_sup);
[x_rel, y_rel] = find(err_rel == min_rel);
% Optima
i0ref_f_l1 = Xa(x_l1,y_l1);
i0ref_a_l1 = Xc(x_l1,y_l1);
i0ref_f_l2 = Xa(x_l2,y_l2);
i0ref_a_l2 = Xc(x_l2,y_l2);
i0ref_f_sup = Xa(x_sup,y_sup);
i0ref_a_sup = Xc(x_sup,y_sup);
i0ref_f_rel = Xa(x_rel, y_rel);
i0ref_a_rel = Xc(x_rel, y_rel);
% Print
fprintf("L1 minimum: i_0,ref,f = %d, i_0,ref,a = %d\n", i0ref_f_l1, i0ref_a_l1);
fprintf("L2 minimum: i_0,ref,f = %d, i_0,ref,a = %d\n", i0ref_f_l2, i0ref_a_l2);
fprintf("Sup minimum: i_0,ref,f = %d, i_0,ref,a = %d\n", i0ref_f_sup, i0ref_a_sup);
fprintf("Rel minimum: i_0.ref.f = %d, i_0,ref,a = %d\n", i0ref_f_rel, i0ref_a_rel);

%% Plot power curves
tab_exp = readtable("Test_850_H2puro_d.csv");
figure;
tcl = tiledlayout(1,2);
% Polarization
nexttile;
plot(Iq, V_obj, 'k');
hold on
plot(Iq, spline(Y_I(:,x_l1,y_l1), V_exp, Iq), '-r');
plot(Iq, spline(Y_I(:,x_l2,y_l2), V_exp, Iq), '-g');
plot(Iq, spline(Y_I(:,x_sup,y_sup), V_exp, Iq), '-b');
plot(Iq, spline(Y_I(:,x_rel,y_rel), V_exp, Iq), '-c');
xlim([0, 2.5]);
ylim([0, 1.3]);
xlabel("Average current [A/cm\textsuperscript{2}]");
ylabel("Voltage [V]");
legend("Exp", "L1", "L2", "Sup", "Rel");
title("Polarization");
% Power
nexttile;
plot(Iq, P_obj, '-k');
hold on
plot(Iq, spline(Y_I(:,x_l1,y_l1), Y_P(:,x_l1,y_l1), Iq), '-r');
plot(Iq, spline(Y_I(:,x_l2,y_l2), Y_P(:,x_l2,y_l2), Iq), '-g');
plot(Iq, spline(Y_I(:,x_sup,y_sup), Y_P(:,x_sup,y_sup), Iq), '-b');
plot(Iq, spline(Y_I(:,x_rel,y_rel), Y_P(:,x_rel,y_rel), Iq), '-c');
xlim([0, 2.5]);
ylim([0, 0.7]);
xlabel("Average current [A/cm\textsuperscript{2}]");
ylabel("Average power [W/cm\textsuperscript{2}]");
legend("Exp", "L1", "L2", "Sup", "Rel");
title("Power");

%% Plot errors
figure;
tcl = tiledlayout(2,2);
% L1
nexttile;
surf(Xa, Xc, err_l1);
xlabel("i_{0,ref,f} [A/cm\textsuperscript{2}]");
ylabel("I_{0,ref,a} [A/cm\textsuperscript{2}]");
title("L^1");
% L2
nexttile;
surf(Xa, Xc, err_l2);
xlabel("I_{0,ref,f} [A/cm\textsuperscript{2}]");
ylabel("I_{0,ref,a} [A/cm\textsuperscript{2}]");
title("L\textsuperscript{2}");
% Sup
nexttile;
surf(Xa, Xc, err_sup);
xlabel("I_{0,ref,f} [A/cm\textsuperscript{2}]");
ylabel("I_{0,ref,a} [A/cm\textsuperscript{2}]");
title("Supremum");
% Rel
nexttile;
surf(Xa, Xc, err_rel);
xlabel("I_{0,ref,f} [A/cm\textsuperscript{2}]");
ylabel("I_{0,ref,a} [A/cm\textsuperscript{2}]");
title("Relative");
% Figure title
title(tcl, "Errors");

%% Plot and export relative error
fig = figure(Position=[300,100,600,600]);
s = surf(Xa, Xc, err_rel);
%set(s, 'FaceColor', [0.4,0.7,0.4]);
clim([0,0.4]);
xlb = xlabel("$i_{0,\mathrm{ref},f}$ [A/cm\textsuperscript{2}]");
ylb = ylabel("$i_{0,\mathrm{ref},a}$ [A/cm\textsuperscript{2}]");
xlb.Position(1) = xlb.Position(1) - 0.15;
xlb.Position(2) = xlb.Position(2) - 0.05;
ylb.Position(1) = ylb.Position(1) + 0.65 ;
ylb.Position(2) = ylb.Position(2) + 2;
title("Relative error");
view(40,35);
% Save plot
ax = gca;
%exportgraphics(ax, "i0ref error.pdf", ContentType="vector");
exportgraphics(ax, "i0ref error.png");