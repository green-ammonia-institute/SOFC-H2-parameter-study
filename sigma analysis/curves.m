% Polarization and power plots

% Models
kinetics = ["BV"];%["BV", "CJ"];
mass = ["MS", "F"];
% Colors
left_color = [.5 .5 0];
right_color = [0 .5 .5];
% Set default plot settings
fsize = 22;
set(0,"defaultTextInterpreter", "latex");
set(0,"defaultLegendInterpreter", "latex");
set(0,'defaultAxesFontSize', fsize);
set(0,'defaultAxesColorOrder',[left_color; right_color]);
set(0,'defaultaxeslinestyleorder',{'-o', '-pentagram', '-square','-v', '-^'}) %or whatever you want
% Axis limits
Ilim = [0, 2.6];
Vlim = [0, 1.6];
Plim = [0, 0.9];

% Load experimental data
tab_exp = readtable("Test_850_H2puro_d.csv");
V_exp = tab_exp{:,1};
I_exp = tab_exp{:,3};
P_exp = tab_exp{:,5};
% Iterate through models
for i = 1:length(kinetics)
    for j = 1:length(mass)
        % Create figure and set colors
        fig = figure(Position=[400,200,1000,500]);
        % Define model
        model = kinetics(i) + "-" + mass(j);
        % Load current data
        tab_current = readtable(model + " - Current.csv");
        % Pre-process data
        [s_ion, V, I] = preprocess_sigma(tab_current);
        % Plot curves for some s_ion
        for k = 1:3
            l = 2*k-1;
            hold on
            xlabel("Average current [A/cm\textsuperscript{2}]");
            xlim(Ilim);
            % Polarization curve
            yyaxis left
            ylabel("Voltage [V]");
            ylim(Vlim);
            plot(I(l,:), V, color=left_color, DisplayName="Voltage, $\sigma_{\mathrm{ion}}$ = "+s_ion(l)+" [S/cm]");
            % Power curve
            yyaxis right
            ylabel("Average power [W/cm\textsuperscript{2}]");
            ylim(Plim);
            plot(I(l,:), I(l,:).*V', color=right_color, DisplayName="Power, $\sigma_{\mathrm{ion}}$ = "+s_ion(l)+" [S/cm]");
        end
        % Plot experimental curves
        yyaxis left
        plot(I_exp, V_exp, Marker="+", color="black", DisplayName="Voltage, experiment");
        yyaxis right
        plot(I_exp, I_exp.*V_exp, Marker="+", color="red", DisplayName="Power, experiment");
        % Legend
        title("Polarization \& power curves")
        legend(FontSize=12);
        % Save plot
        ax = gca;
        exportgraphics(ax, model + " - Polarization.pdf", BackgroundColor="none", ContentType="vector");
        exportgraphics(ax, model + " - Polarization.png", BackgroundColor="none");
        % Calculate maximum power variation
        maxP = max(I.*V', [], 2);
        message = "The max. power output of model %s at s_ion = %d is %f\n";
        fprintf(message, model, s_ion(1), maxP(1));
        fprintf(message, model, s_ion(end), maxP(end));
        fprintf("Max. power output increased by %f percent\n", 100*abs(maxP(end)-maxP(1))/maxP(1));
    end
end

%% Plot power vs sigma_ion
% Models
kinetics = ["BV"];%["BV", "CJ"];
mass = ["MS", "F"];
% Set default plot settings
fsize = 12;
set(0,'defaultAxesFontSize', fsize);

% Iterate through models
for i = 1:length(kinetics)
    for j = 1:length(mass)
        % Define model
        model = kinetics(i) + "-" + mass(j);
        % Load current data
        tab_current = readtable(model + " - Current.csv");
        % Pre-process data
        [s_ion, V, I] = preprocess_sigma(tab_current);
        nV = length(V);
        % Plot
        fig = figure;
        colors = jet(nV);
        %ylim(Plim);
        xlabel("Ionic conductivity [S/cm]");
        ylabel("Average power [W/cm\textsuperscript{2}]");
        for k = 1:nV
            xlabel("Ionic conductivity [S/cm]");
            ylabel("Average power [W/cm\textsuperscript{2}]");
            semilogx(s_ion, I(:,k) * V(k), Marker="*", Color=colors(k,:),DisplayName="$V_{\mathrm{cell}}$ = "+V(k)+" [V]");
            hold on
        end
        % Legend
        title("Power vs. $\sigma_{\mathrm{ion}}$")
        legend(FontSize=10, Location="eastoutside");
        ax = gca;
        exportgraphics(ax, model + " - Power vs sigma.pdf", BackgroundColor="none", ContentType="vector");
    end
end

%% Convergence
% Error plots

% Models
kinetics = ["BV"];%, "CJ"];
mass = ["MS", "F"];
% Set default plot settings
fsize = 23;
set(0,'defaultAxesFontSize', fsize);
% Colors
set(0,'defaultAxesColorOrder', [[1 0.9 0]; [0.5 0 0.5]; [0.2 1 0.7]; [0.5 0.5 1]; [0 1 0]]);
%set(0,'defaultaxeslinestyleorder',{'-o', '-pentagram', '-square','-v', '-^'}) %or whatever you want
% Axis limits
names = {"Current", "Momentum - Air", "Mass - Air", "Momentum - Fuel", "Mass - Fuel"};

% Iterate through models
for i = 1:length(kinetics)
    for j = 1:length(mass)
        % Create figure and set colors
        fig = figure(Position=[100,100,1400,600]);
        % Define model
        model = kinetics(i) + "-" + mass(j);
        % Load current data
        tab_error = readtable(model + " - Convergence.csv");
        % Plot curves
        for k = 1:5
            xlabel("Iteration");
            ylabel("Error");
            semilogy(tab_error{:,1}, tab_error{:,k+1}, '-', DisplayName=names{k});
            hold on
        end
        ax = gca;
        %ax.FontSize = 15;
        % Legend
        title("Segregated solver error")
        legend(FontSize=14, Location="southwest");
        % Save plot
        exportgraphics(ax, model + " - Convergence.pdf", BackgroundColor="none", ContentType="vector");
    end
end

%% Log
% Polarization and curve plots

% Models
kinetics = ["BV", "CJ"];
mass = ["MS", "F"];
% Colors
left_color = [.5 .5 0];
right_color = [0 .5 .5];
set(0,'defaultAxesColorOrder',[left_color; right_color]);
set(0,'defaultaxeslinestyleorder',{'-o', '-pentagram', '-square','-v', '-^'}) %or whatever you want
% Axis limits
Ilim = [-6, 1];
Vlim = [0, 1.3];
Plim = [0, 0.61];

% Load experimental data
tab_exp = readtable("Test_850_H2puro_d.csv");
V_exp = tab_exp{:,1};
I_exp = tab_exp{:,3};
P_exp = tab_exp{:,5};
logI_exp = log(I_exp);
% Iterate through models
for i = 1:length(kinetics)
    for j = 1:length(mass)
        % Create figure and set colors
        fig = figure(Position=[400,200,1000,500]);
        % Define model
        model = kinetics(i) + "-" + mass(j);
        % Load current data
        tab_current = readtable(model + " - Current.csv");
        % Pre-process data
        [s_ion, V, I] = preprocess_data(tab_current);
        % Plot curves for some s_ion
        for k = 1:3
            l = 2*k-1;
            hold on
            logI = log(I(l,:));
            P = I(l,:).*V';
            xlabel("Average current [A/cm\textsuperscript{2}]");
            xlim(Ilim);
            % Polarization curve
            yyaxis left
            ylabel("Voltage [V]");
            ylim(Vlim);
            plot(logI, V, color=left_color, DisplayName="Voltage, $\sigma_{\mathrm{ion}}$ = "+s_ion(l));
            % Power curve
            yyaxis right
            ylabel("Average power [W/cm\textsuperscript{2}]");
            ylim(Plim);
            plot(logI, P, color=right_color, DisplayName="Power, $\sigma_{\mathrm{ion}}$ = "+s_ion(l));
            % Terminal slope
            slope = (V(end) - V(end-1))/(logI(end) - logI(end-1)); 
            fprintf("sigma_ion = %d, slope = %d\n", s_ion(l), slope);
        end
        % Plot experimental curves
        yyaxis left
        plot(log(I_exp), V_exp, Marker="+", color="black", DisplayName="Voltage, experiment");
        yyaxis right
        plot(log(I_exp), I_exp.*V_exp, Marker="+", color="red", DisplayName="Power, experiment");
        % Terminal slope
        slope = (V_exp(end) - V_exp(end-1))/(logI_exp(end) - logI_exp(end-1));
        fprintf("Experimental slope = %d\n", s_ion(l), slope);
        % Legend
        title("Polarization & power curves")
        legend(FontSize=8, Location="southwest");
        % Save plot
        ax = gca;
        exportgraphics(ax, model + " - LogPolarization.pdf", BackgroundColor="none");
    end
end

%% Mole fraction plots
%Plot settings
colormap("jet");
% Load data
tab_xH2 = readtable("BV-MS - xH2.txt");
X = tab_xH2{:,1};
Z = tab_xH2{:,2};
xH2 = tab_xH2{:,3};
% Interpolate
F = scatteredInterpolant(X, Z, xH2);
% Create grid and plot
a = linspace(min(X), max(X), 100);
b = linspace(min(Z), max(Z), 70);
[XX,ZZ] = meshgrid(a,b)
P = [XX(:) ZZ(:)];
inter = F(P);
inter = reshape(inter, 100, 70);
%s = contourf(XX, ZZ, inter, 1000, LineColor='none');
s = imagesc(a, b, inter);
clim([0,1])
colorbar;