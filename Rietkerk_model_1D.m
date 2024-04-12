% 1D_Rietkerk_model
tic 
% clear
% Constants:
g = 0.05; 
c = 10;
k1 = 5;
k2 = 5;
d = 0.25;
alpha = 0.2;
rW = 0.2;
w0 = 0.2;
DP = 0.1;
DW = 0.1;
DO = 100; 

% P biomass, W Soil water, O Surface water

% Vary rainfall:
R = 1.2;

% time step
len = 100;
dx = 1;
N = len/dx;
T = 4000; % total time 
dt = 0.0005; % time step
n = floor(T / dt); % number of iterations
 % 1:n/8:n


% Initial Conditions
% start near an equilibria

% For R >= 1 
% Uniformly vegetated state
W = d*k1/(c*g - d);
P = (R-(rW.*W)).*(W+k1)./(g*W);
O = R/alpha * (P + k2)./(P + (k2*w0));


xplot = 0:dx:len-dx;


L = 4*2*pi/len; % Of the form n*2*pi/len
sigma = h11(R, L);



pert = ones(3,1) .* real(exp(1i*sigma.*dt).*sin(L*xplot)); 


PWO = [P;W;O] + (pert.* 1);

 % PWO = PWO_IC; % should be at a non-homogeneous steady state



% Combined laplcian
lap_cbn_1d = @(Z)  (Z(1:3,1:end-2) + Z(1:3,3:end) - 2 * Z(1:3,2:end-1)) / dx^2 ;


% Combined homogeneous Rietkerk model
R2002 = @(PWO) [c*g*PWO(2,:).*PWO(1,:)./(PWO(2,:)+k1) - d*PWO(1,:); ...
                alpha*PWO(3,:).*(PWO(1,:)+(k2*w0))./(PWO(1,:)+k2) - g*PWO(2,:).*PWO(1,:)./(PWO(2,:)+k1) - rW*PWO(2,:); ...
                R - alpha*PWO(3,:).*(PWO(1,:)+(k2*w0))./(PWO(1,:)+k2)];


j = 0; % For subplot loop

figure
for i = 1:n

    % Laplacian
    deltaPWO = lap_cbn_1d(PWO);
    
    % inside the grid
    PWOc = [PWO(1, 2:end-1); PWO(2, 2:end-1); PWO(3, 2:end-1)];
    
    % Equations
    PWO(1:3, 2:end-1) = PWOc + dt.*R2002(PWOc) + dt.*[DP;DW;DO].*deltaPWO;

    % % Neumann BCs:
    % PWO(:,1) = PWO(:,2);
    % PWO(:,end) = PWO(:,end-1);

    % Periodic BC's:
    PWO(:,1) = PWO(:,end-1);
    PWO(:,end) = PWO(:,2); 

    if ismember(i, [1,100/dt, 400/dt, 4000/dt])
        j = j+1;
        subplot(1,4, j)
        plot(xplot,PWO(1,:), 'r', 'LineWidth', 3);
        hold on
        plot(xplot,PWO(2,:), 'g', 'LineWidth', 3);
        plot(xplot,PWO(3,:), 'b', 'LineWidth', 3);
        title(['t = ', num2str(i * dt)]);
        set(gca, 'FontSize', 15)
        xlabel("Location [m]")
        % ylabel("Plant density [g/m^2], surface/soil water [mm]")
        % legend('Plant biomass, P', 'Soil water, W', 'Surface water, O')
        ylim([-0.1, 20])

    end
end

subplot(1,4,1)
ylabel({'Plant density [g/m^2],'; 'surface/soil water [mm]'}, "FontSize", 15)

%%
% Large final plot at time T
figure
plot(xplot,PWO(1,:), 'r', 'LineWidth', 3);
hold on
plot(xplot,PWO(2,:), 'g', 'LineWidth', 3);
plot(xplot,PWO(3,:), 'b', 'LineWidth', 3);
title(['t = ', num2str(i * dt)]);
subtitle(['R = ', num2str(R)]);
set(gca, 'FontSize', 15)
xlabel("Spatial localtion [m]")
ylabel("Plant density [g/m^2], surface/soil water [mm]")
legend('Plant biomass, P', 'Soil water, W', 'Surface water, O')
ylim([-0.1, 20])

PWO_IC = PWO; % store the last values


toc
% Equilibria plots
syms Rain

Wneq = d*k1/(c*g - d);
Pneq = @(Rain) (Rain-(rW*Wneq))*(Wneq+k1)/(g*Wneq);
Oneq = @(Rain) Rain/alpha * (Pneq(Rain) + k2)/(Pneq(Rain) + (k2*w0));

Pzero = 0;
Wzero = @(Rain) Rain/(rW);
Ozero = @(Rain) Rain/(alpha*w0);

Biomass = trapz(PWO(1,:));
Avg_Biomass = Biomass/len;

% Bifurcation diagram for biomass
figure
fplot(@(Rain)Pneq(Rain), [0, 3], 'LineWidth', 2, 'Color', 'b')
hold on
fplot(@(Rain)Pneq(Rain), [k1*rW*d/(c*g - d), max_unstab_R], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'w')
fplot(Pzero, [0,3], 'LineWidth', 2, 'Color', 'b')
fplot(Pzero, [k1*rW*d/(c*g - d), 3], 'LineWidth', 2, 'Color', 'w')
fplot(Pzero, [k1*rW*d/(c*g - d), 3], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r')

plot(R, max(PWO(1,:)), 'x', 'MarkerSize', 20, 'LineWidth', 3, 'Color', 'r')
plot(R, Avg_Biomass, 'x', 'MarkerSize', 20, 'LineWidth', 3, 'Color', 'g')
xline(0)
xlim([0.5, 1.5])
ylim([0,20])
set(gca, 'FontSize', 15)
ylabel('Biomass, P', 'FontSize', 24)
xlabel('Rainfall, R', 'FontSize', 24)

%% test
plot(0.8, max_P_tempsave, 'x', 'MarkerSize', 20, 'LineWidth', 3, 'Color', 'r')
plot(0.8, Biomass_tempsave, 'x', 'MarkerSize', 20, 'LineWidth', 3, 'Color', 'g')





%% Equilibria plots
% Soil water
figure
fplot(Wneq, [k1*rW*d/(c*g - d), 3], 'LineWidth', 2, 'Color', 'b')
hold on
fplot(@(Rain)Wzero(Rain), [k1*rW*d/(c*g - d), 3], 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--')
fplot(@(Rain)Wzero(Rain), [0, k1*rW*d/(c*g - d)], 'LineWidth', 2, 'Color', 'r')
set(gca, 'FontSize', 15)
ylabel('Soil water, W [mm]', 'FontSize', 24)
xlabel('Rainfall, R [mm/day]', 'FontSize', 24)
ylim([0,10])
% legend('P \neq 0', 'P = 0')

% Surface water
figure 
fplot(@(Rain)Oneq(Rain), [k1*rW*d/(c*g - d), 3], 'LineWidth', 2, 'Color', 'b')
hold on
fplot(@(Rain)Ozero(Rain), [0, k1*rW*d/(c*g - d)], 'LineWidth', 2, 'Color', 'r')
fplot(@(Rain)Ozero(Rain), [k1*rW*d/(c*g - d), 3], 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--')
set(gca, 'FontSize', 15)
ylabel('Surface water, O [mm]', 'FontSize', 24)
xlabel('Rainfall, R [mm/day]', 'FontSize', 24)
% legend('P \neq 0', 'P = 0')
ylim([0,50])

% Plant Biomass
figure
fplot(@(Rain)Pneq(Rain), [0, 3], 'LineWidth', 2, 'Color', 'b')
hold on
fplot(Pzero, [0,k1*rW*d/(c*g - d)], 'LineWidth', 2, 'Color', 'r')
fplot(Pzero, [k1*rW*d/(c*g - d), 3], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r')
set(gca, 'FontSize', 15)
ylabel('Biomass, P [g/m^2]', 'FontSize', 24)
xlabel('Rainfall, R [mm/day]', 'FontSize', 24)
ylim([0,50])

%% Water infiltration plot
alpha = 0.2;
O = 5; % arb
k2 = 5;
w0 = 0.2;

% Define Inf as an anonymous function
Inf = @(Plant) alpha * O * (Plant + k2 * w0) ./ (Plant + k2);

% Plot Inf
figure
fplot(@(Plant) Inf(Plant), [0, 200], 'Color', 'k', 'LineWidth', 2)
hold on
set(gca, 'FontSize', 15);

% Asymptote
asymp = w0 * O; 
yline(asymp, 'LineStyle', ':', 'Color', 'b', 'LineWidth', 2);

ylim([0,asymp*1.1])
xticks([])
yticks([alpha * w0 * O, w0 * O])
yticklabels({'\alphaOW_0', '\alphaO'})

xlabel('Plant Density', 'FontSize', 24)
ylabel('Water Infiltration', 'FontSize', 24)


%% Water uptake plot
g = 0.05;
k1 = 5;
P = 20;
c = 10;

uptake = @(Soil) P * g * Soil ./ (Soil+k1);


figure
fplot(@(Soil) uptake(Soil), [0, 200], 'Color', 'k', 'LineWidth', 2)
hold on

% Asymptote
asymp = P*g; 
yline(asymp, 'LineStyle', ':', 'Color', 'b', 'LineWidth', 2);
ylim([0,asymp*1.1])

xticks([])
yticks(P*g)

yticklabels({'Pg or Pcg'})
set(gca, 'FontSize', 15)

xlabel('Soil Water Availability', 'FontSize', 24)
ylabel('Plant Water Uptake', 'FontSize', 24)




