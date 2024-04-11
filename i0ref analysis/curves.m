% Polarization and power plots

% Models
kinetics = ["BV"];
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
        [i0refs, V, I] = preprocess_current(tab_current);
        % Plot curves for i_0,ref,f
        for k = 1:length(i0refs)
            hold on
            xlabel("Average current [A/cm\textsuperscript{2}]");
            xlim(Ilim);
            % Polarization curve
            yyaxis left
            ylabel("Voltage [V]");
            ylim(Vlim);
            plot(I(k,:), V, color=left_color, DisplayName="Voltage, $i_{0,\mathrm{ref},f}$ = "+i0refs(k)+" [A/cm\textsuperscript{2}]");
            % Power curve
            yyaxis right
            ylabel("Average power [W/cm\textsuperscript{2}]");
            ylim(Plim);
            plot(I(k,:), I(k,:).*V', color=right_color, DisplayName="Power, $i_{0,\mathrm{ref},f}$ = "+i0refs(k)+" [A/cm\textsuperscript{2}]");
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
        message = "The max. power output of model %s at i0ref_f = %d is %f\n";
        fprintf(message, model, i0refs(1), maxP(1));
        fprintf(message, model, i0refs(end), maxP(end));
        fprintf("Max. power output increased by %f percent\n", 100*abs(maxP(end)-maxP(1))/maxP(1));
    end
end

%% Plot power vs sigma_ion
% Models
kinetics = ["BV"];
mass = ["MS", "F"];
% Set default plot settings
fsize = 12;
set(0,"defaultTextInterpreter", "latex");
set(0,"defaultLegendInterpreter", "latex");
set(0,'defaultAxesFontSize', fsize);

% Iterate through models
for i = 1:length(kinetics)
    for j = 1:length(mass)
        % Define model
        model = kinetics(i) + "-" + mass(j);
        % Load current data
        tab_current = readtable(model + " - Current.csv");
        % Pre-process data
        [i0refs, V, I] = preprocess_current(tab_current);
        nV = length(V);
        % Plot
        fig = figure;
        colors = jet(nV);
        %ylim(Plim);
        %ylabel("Average power [W/cm^2]");
        for k = 1:nV
            xlabel("Reference exchange current density [A/cm\textsuperscript{2}]");
            ylabel("Average power [W/cm\textsuperscript{2}]");
            semilogx(i0refs, I(:,k) * V(k), Marker="*", Color=colors(k,:),DisplayName="$V_{\mathrm{cell}}$ = "+V(k)+" [V]");
            hold on
        end
        % Legend
        title("Power vs. $i_{0,\mathrm{ref},f}$")
        legend(FontSize=10, Location="eastoutside");
        ax = gca;
        exportgraphics(ax, model + " - Power vs i0ref.pdf", BackgroundColor="none", ContentType="vector");
    end
end

%% Convergence
% Error plots

% Models
kinetics = ["BV"];
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
        % Change tick font size
        ax = gca;
        %ax.FontSize = 15;
        % Legend
        title("Segregated solver error")
        legend(FontSize=14, Location="southwest");
        % Save plot
        exportgraphics(ax, model + " - Convergence.pdf", BackgroundColor="none", ContentType="vector");
    end
end