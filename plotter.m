% Get current working folder info
currentFolder = pwd;
currentFolder = strcat(currentFolder, '/CDC_2023_Figures/');
set(groot,'defaultLegendInterpreter','latex');

% Set the time horizon
T = 30;
tvec = 0:T-1;

% Load all the data
load('AdverseDistHinftyData.mat');
load('AdverseDistMinimaxData.mat');
load('ConfuseDistHinftyData.mat');
load('ConfuseDistMinimaxData.mat');
load('SineDistHinftyData.mat');
load('SineDistMinimaxData.mat');
load('RegretAdverseData.mat');
load('RegretConfuseData.mat');
load('RegretSineData.mat');

%% Plot the control inputs in all three cases
figNum = 1;
figure(figNum);
set(gcf, 'Position', get(0, 'Screensize'));
hold on;
grid on;
plot(tvec, sine_dist_u_minmax(1,:)' - sine_dist_u_hinfty(1,:)', '-*');
plot(tvec, adverse_dist_u_minmax(1,:)' - adverse_dist_u_hinfty(1,:)', '-*');
plot(tvec, confuse_dist_u_minmax(1,:)' - confuse_dist_u_hinfty(1,:)', '-*');
xlabel('time');
ylabel('controls');
xlim([0, T]);
ylim([-26, 12]);
legend('$u^{\dagger}_{\mathrm{sine}} - u^{\star}_{\mathrm{sine}}$', '$u^{\dagger}_{\mathrm{adv}} - u^{\star}_{\mathrm{adv}}$', '$u^{\dagger}_{\mathrm{thm} 1} - u^{\star}_{\mathrm{thm} 1}$', 'Location','southeast');
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 6);
set(a, 'linewidth', 6);
set(a, 'FontSize', 60);
saveas(figure(figNum), [currentFolder 'Controls.png']);

%% Plot the difference of states in all three cases
figNum = figNum + 1;
tvec = 0:T;
figure(figNum);
set(gcf, 'Position', get(0, 'Screensize'));
hold on;
grid on;
plot(tvec, sine_dist_x_minmax(1,:)' - sine_dist_x_hinfty(1,:)', '-*');
plot(tvec, adverse_dist_x_minmax(1,:)' - adverse_dist_x_hinfty(1,:)', '-*');
plot(tvec, confuse_dist_x_minmax(1,:)' - confuse_dist_x_hinfty(1,:)', '-*');
xlabel('time');
ylabel('states');
xlim([0, T]);
ylim([-15, 10]);
legend('$x^{\dagger}_{\mathrm{sine}} - x^{\star}_{\mathrm{sine}}$', '$x^{\dagger}_{\mathrm{adv}} - x^{\star}_{\mathrm{adv}}$', '$x^{\dagger}_{\mathrm{thm} 1} - x^{\star}_{\mathrm{thm} 1}$', 'Location','southeast');
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 6);
set(a, 'linewidth', 6);
set(a, 'FontSize', 60);
saveas(figure(figNum), [currentFolder 'States.png']);

%% Plot the regrets in all three cases
figNum = figNum + 1;
tvec = 0:T;
figure(figNum);
set(gcf, 'Position', get(0, 'Screensize'));
hold on;
grid on;
plot(tvec, instSineRegret, '-*');
plot(tvec, cumSineRegret, '-*');
plot(tvec, instAdvRegret, '-*');
plot(tvec, cumAdvRegret, '-*');
plot(tvec, instConfRegret, '-*');
plot(tvec, cumConfRegret, '-*');
xlabel('time');
ylabel('regret');
xlim([0, T]);
ylim([0, 3000]);
legend('$\mathcal{R}_{\mathrm{sine}}$', '$\tilde{\mathcal{R}}_{\mathrm{sine}}$', '$\mathcal{R}_{\mathrm{adv}}$', '$\tilde{\mathcal{R}}_{\mathrm{adv}}$', '$\mathcal{R}_{\mathrm{thm} 1}$', '$\tilde{\mathcal{R}}_{\mathrm{thm} 1}$','Location','northeast', 'NumColumns',3);
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 6);
set(a, 'linewidth', 6);
set(a, 'FontSize', 60);
saveas(figure(figNum), [currentFolder 'Regret.png']);