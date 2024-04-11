% Cut line plots

% Lengths [um]
H_gde = 50;
H_electrolyte = 200;
H_tube = 2000;
H_tot = 2*H_tube + H_electrolyte;
% Reference exchange currents [A/cm^2]
sigmas = [1e-3, 5e-3, 1e-2, 5e-2, 1e-1];
ns = length(sigmas);
% Models
kinetics = ["BV"];%["BV", "CJ"];
mass = ["MS", "F"];
% Axis limits
Llim = [0, H_tot];
MSH2lim = [0.79, 0.83];
MSO2lim = [0.175, 0.21];
FLH2lim = [0.996, 1];
FLO2lim = [0.2075, 0.21];
ylims = zeros(2,2,2);
ylims(1,1,:) = MSH2lim;
ylims(1,2,:) = MSO2lim;
ylims(2,1,:) = FLH2lim;
ylims(2,2,:) = FLO2lim;
% Set default plot settings
colors = [0 0 1; 0 1 0; 1 0 0];
fsize = 18;
set(0,'defaultAxesFontSize', fsize);
set(0,"defaultTextInterpreter", "latex");
set(0,"defaultLegendInterpreter", "latex");

% Iterate through models
for i = 1:length(kinetics)
    for j = 1:length(mass)
        % Create figure and set colors
        fig = figure(Position=[400,200,1000,500]);
        % Define model
        model = kinetics(i) + "-" + mass(j);
        % Load concentration data
        tab_H2 = readtable(model + " - Cut Line H2.csv");
        tab_O2 = readtable(model + " - Cut Line O2.csv");
        L = tab_H2{:,1}*1e6;
        xH2 = tab_H2{:,2:end};
        xO2 = tab_O2{:,2:end};
        % Plot curves for i_0,ref,f
        tcl = tiledlayout(2,1, TileSpacing="tight", Padding="compact");
        ax_H2 = nexttile(1);
        ax_O2 = nexttile(2);
        hold(ax_H2);
        hold(ax_O2);
        box off;
        for k = 1:3
            box off;
            m = 2*k-1;
            % Plot H2 curves
            plot(L, xH2(:,m), Parent=ax_H2, Marker="none", Color=colors(k,:), DisplayName="H$_2$, $\sigma_{\mathrm{ion}}$ = "+sigmas(m)+" [S/cm]");
            legend(ax_H2, Location="northwest", FontSize=10);
            set(ax_H2, 'XTick', []);
            % Plot O2 curves
            plot(L, xO2(:,m), Parent=ax_O2, Marker="none", Color=colors(k,:), DisplayName="O$_2$, $\sigma_{\mathrm{ion}}$ = "+sigmas(m)+" [S/cm]");
            legend(ax_O2, Location="northeast", FontSize=10);
        end
        linkaxes([ax_H2, ax_O2],'x');
        % Plot guiding lines
        for l = 1:2
            xline(nexttile(l), H_tube-H_gde, LineStyle="--", Alpha=0.3, HandleVisibility="off");
            xline(nexttile(l), H_tube, LineStyle="--", Alpha=0.3, HandleVisibility="off");
            xline(nexttile(l), H_tube+H_electrolyte, LineStyle="--", Alpha=0.3, HandleVisibility="off");
            xline(nexttile(l), H_tube+H_electrolyte+H_gde, LineStyle="--", Alpha=0.3, HandleVisibility="off");
        end
        % Legend
        xlabel(tcl, "Arc length [$\mu$m]", FontSize=fsize, Interpreter="latex");
        xlim(ax_O2,Llim);
        ylim(ax_H2, ylims(j,1,:));
        ylim(ax_O2, ylims(j,2,:));
        ylabel(tcl,"Mole fraction", FontSize=fsize, Interpreter="latex");
        title(tcl, "Molar fractions, central cut line", FontSize=fsize, Interpreter="latex");
        % Save plot
        exportgraphics(fig, model + " - Cut Line.pdf", BackgroundColor="none", ContentType="vector");
    end
end

%% Old version using stackedplot, ignore

% Iterate through models
for i = 1:length(kinetics)
    for j = 1:length(mass)
        % Create figure and set colors
        fig = figure(Position=[400,200,1000,500]);
        % Define model
        model = kinetics(i) + "-" + mass(j);
        % Load concentration data
        tab_H2 = sortrows(readtable(model + " - Cut Line H2.csv"));
        tab_O2 = sortrows(readtable(model + " - Cut Line O2.csv"));
        L = tab_H2{:,1}*1e6;
        xH2 = tab_H2{:,2:end};
        xO2 = tab_O2{:,2:end};
        tab_x = table(L, xH2, xO2);
        % Plot curves for i_0,ref,f
        s = stackedplot(tab_x, "XVariable","L");
        axesProps = struct(s.AxesProperties(2));  % using struct() is undocumented
        axesProps.Axes.XLabel.Interpreter = 'tex';
        s.XLimits = Llim;
        s.XLabel = "Arc length [\mu m]";
        s.Title = "Molar fractions, central cut line";
        % Set H2 plot properties
        %axesProps.LegendLabels_I.Interpreter = 'tex';
        s.AxesProperties(1).LegendLabels = arrayfun(@(i) "H_2, i_{0,ref,f} = "+i+" [A/cm^2]", i0refs);%["H_2, i_{0,ref,f} = ", "b", "c"];
        s.AxesProperties(1).YLimits = H2lim;
        % Set O2 plot properties
        % Legend
        %legend(FontSize=8);
        % Save plot
        %ax = gca;
        %exportgraphics(ax, model + " - Polarization.pdf", BackgroundColor="none");
    end
end
