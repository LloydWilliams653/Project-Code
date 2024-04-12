% 1D_Rietkerk_model plotting the non-homogeneous solutions
% Need to run stability analysis code prior to this
% len = 100 & Rainfall_drop = 0.01 ~ 13.5 hours
% len = 100 & Rainfall_drop = 0.40 ~ 25 minutes
% len = 200 & Rainfall_drop = 0.01 ~ 17 hours
% len = 200 & Rainfall_drop = 0.40 ~ 33 minutes


tic 

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

% Rainfall:
Initial_rain = 1.2;

R = Initial_rain;


% time step
len = 200; % Domain size [m]
dx = 1;
N = len/dx;
T = 4000; % total time 
dt = 0.0005; % time step
n = floor(T / dt); % number of iterations


% Initial Conditions
% start near an equilibria
% For R >= k1*rW*d/(c*g - d) (>= 1)
% Uniformly vegetated state
W = d*k1/(c*g - d);
P = (R-(rW.*W)).*(W+k1)./(g*W);
O = R/alpha * (P + k2)./(P + (k2*w0));


xplot = 0:dx:len-dx;


% 0.01, 0.4
Rainfall_drop = 0.4; % Change in rainfall after T

% Range of R value to vary
min_R = 0;
max_R = 1.4;
samples = 1000;% number of samples of biomass taken at each rainfall


% Decreasing R
J = ceil(abs(R-min_R)/Rainfall_drop);
 % J = 4; % test
 
% Increasing R from 1.25 to 1.4
J2 = ceil(abs(R-max_R)/Rainfall_drop);
% J2 = 4; % test


% Spatial frequency n = 2;
L =  2*2*pi/len;
sigma = h11(R, L);


% Zero vectors for the outputs of Max plant peaks, average biomass & Biomass distribution 
max_P = zeros(2,J+J2-2);
Biomass = zeros(2,J+J2-2);
m = 0; % For store_bio loop
store_bio = zeros(J*(samples),N);



% Spatially inhomogeneous Perturbation
pert = ones(3,1) .* real(exp(1i*sigma.*dt).*sin(L*xplot)); 


% loop for increasing R
for j2 = 1:J2+1

% Set the initial point
    if j2 == 1
  PWO = [P;W;O] + (pert.* 1);

    else
  PWO = PWO_IC;
    end



% Laplcian equation (Diffusion)
lap_cbn_1d = @(Z)(Z(1:3,1:end-2) + Z(1:3,3:end) - 2*Z(1:3,2:end-1))/dx^2;


% Homogeneous Rietkerk model
R2002 = @(PWO) [c*g*PWO(2,:).*PWO(1,:)./(PWO(2,:)+k1) - d*PWO(1,:); ...
                alpha*PWO(3,:).*(PWO(1,:)+(k2*w0))./(PWO(1,:)+k2) - g*PWO(2,:).*PWO(1,:)./(PWO(2,:)+k1) - rW*PWO(2,:); ...
                R - alpha*PWO(3,:).*(PWO(1,:)+(k2*w0))./(PWO(1,:)+k2)];



% Euler loop
for i = 1:n

    % Laplacian
    deltaPWO = lap_cbn_1d(PWO);
    
    % inside the grid
    PWOc = [PWO(1, 2:end-1); PWO(2, 2:end-1); PWO(3, 2:end-1)];
    
    % Equations
    PWO(1:3, 2:end-1) = PWOc + dt.*R2002(PWOc) + dt.*[DP;DW;DO].*deltaPWO;


    % Periodic BC's:
    PWO(:,1) = PWO(:,end-1);
    PWO(:,end) = PWO(:,2); 

if j2 > 1 % Exclude the uniform to equilibria transition
    if mod(i,n/samples) == 0
        m = m+1;
        store_bio(m,:) = PWO(1,:);
    end
end

end


PWO_IC = PWO; % store the last values to be the new IC


if round(R,2) == 1.2
    PWO_IC_120 = PWO; % Store PWO for R = 1.2
end

% So we omit the perturbed homogeneous transition to non-homog. equilibria
if j2 > 1
max_P(:,j2-1) = [R;max(PWO(1,:))]; % store max peak
Biomass(:,j2-1) = [R;trapz(PWO(1,:))]; % store area under the biomass curve
end

R = round(R + Rainfall_drop,2); % Increase R for the next iteration


end

% invert the max_P
max_P_pt1 = max_P(:, 1:j2-1);
max_P_pt2 = max_P(:, j2+1:end);
invert_max_P_pt1 = max_P_pt1(:, end:-1:1);
max_P = [invert_max_P_pt1, max_P_pt2];


% invert the biomass
Biomass_pt1 = Biomass(:, 1:j2-1);
Biomass_pt2 = Biomass(:, j2+1:end); % rest of the zeroes
invert_Biomass_pt1 = Biomass_pt1(:, end:-1:1);
Biomass = [invert_Biomass_pt1, Biomass_pt2];

% invert the store_bio
store_bio_pt1 = store_bio(1:m, :);
store_bio_pt2 = store_bio(m+1:end, :); % rest of the zeroes
invert_store_bio_pt1 = store_bio_pt1(end:-1:1, :);
store_bio = [invert_store_bio_pt1; store_bio_pt2];


% Reset the new rainfall, back to it's initial value
R = round(Initial_rain + Rainfall_drop,2);




for j = 1:J+1

% Set the initial point
    if j == 1
 PWO = PWO_IC_120;
    else
PWO = PWO_IC;
    end

R = round(R - Rainfall_drop,2); % Reduce R for the next iteration

if R < 0 % Ensure no negative rainfalls
    R = 0;
end


% Laplcian (Diffusion)
lap_cbn_1d = @(Z)  (Z(1:3,1:end-2) + Z(1:3,3:end) - 2 * Z(1:3,2:end-1)) / dx^2 ;


% Homogeneous Rietkerk model
R2002 = @(PWO) [c*g*PWO(2,:).*PWO(1,:)./(PWO(2,:)+k1) - d*PWO(1,:); ...
                alpha*PWO(3,:).*(PWO(1,:)+(k2*w0))./(PWO(1,:)+k2) - g*PWO(2,:).*PWO(1,:)./(PWO(2,:)+k1) - rW*PWO(2,:); ...
                R - alpha*PWO(3,:).*(PWO(1,:)+(k2*w0))./(PWO(1,:)+k2)];



% Euler loop
for i = 1:n

    % Laplacian
    deltaPWO = lap_cbn_1d(PWO);

    % inside the grid
    PWOc = [PWO(1, 2:end-1); PWO(2, 2:end-1); PWO(3, 2:end-1)];

    % Equations
    PWO(1:3, 2:end-1) = PWOc + dt.*R2002(PWOc) + dt.*[DP;DW;DO].*deltaPWO;

    % Periodic BC's:
    PWO(:,1) = PWO(:,end-1);
    PWO(:,end) = PWO(:,2);


    if mod(i,n/samples) == 0
        m = m+1;
        store_bio(m,:) = PWO(1,:); 
    end

if round(R,2) == 0.8
    PWO_IC_08 = PWO;
end

end
PWO_IC = PWO; % store the last values to be the new IC


max_P(:,j2+j-1) = [R;max(PWO(1,:))]; % store max peak
Biomass(:,j2+j-1) = [R;trapz(PWO(1,:))]; % store area under the biomass curve

end

%%
num_colors = 256;
custom_colormap = [linspace(1,0,num_colors)', linspace(1, 0.5, num_colors)', linspace(1,0,num_colors)']; 

% Biomass distribution plot
figure
imagesc(store_bio)
colormap(custom_colormap)
axis xy
colorbar 
set(gca, 'FontSize', 15)
ylabel('Time [day]',  'FontSize', 24)
xlabel('Spatial location [m]',  'FontSize', 24)
ylabel(colorbar, 'Biomass Density [g/m^2]', 'FontSize', 18)

% Change the y axis labels, so it displays time
yticklabs = get(gca, 'yticklabels');
ytickvals = str2double(yticklabs);

ytickvals_new = T/samples * ytickvals;
yticklabels_new = cellstr(num2str(ytickvals_new(:)));

set(gca, 'yticklabels', yticklabels_new);


% Equilibria plots
syms Rain

Wneq = d*k1/(c*g - d);
Pneq = @(Rain) (Rain-(rW*Wneq))*(Wneq+k1)/(g*Wneq);
Oneq = @(Rain) Rain/alpha * (Pneq(Rain) + k2)/(Pneq(Rain) + (k2*w0));

Pzero = 0;
Wzero = @(Rain) Rain/(rW);
Ozero = @(Rain) Rain/(alpha*w0);

% Plotting average biomass, & biomass homogeneous & non-homogeneous equilibria 
figure
fplot(@(Rain)Pneq(Rain), [0, 3], 'LineWidth', 2, 'Color', 'b')
hold on

fplot(@(Rain)Pneq(Rain), [k1*rW*d/(c*g - d), max_unstab_R], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'w') 

fplot(Pzero, [0,3], 'LineWidth', 2, 'Color', 'b')
fplot(Pzero, [k1*rW*d/(c*g - d), 3], 'LineWidth', 2, 'Color', 'w')
fplot(Pzero, [k1*rW*d/(c*g - d), 3], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r')


% Maximum biomass peaks
plot(max_P(1,:), max_P(2,:), '-x', 'MarkerSize', 10, 'LineWidth', 3, 'Color', 'r')
xlim([0.4, 1.5])
ylim([0,20])
set(gca, 'FontSize', 15)
ylabel('Biomass, P [g/m^2]',  'FontSize', 24)
xlabel('Rainfall, R, [mm/day]',  'FontSize', 24)

% Average Biomas
plot(Biomass(1,:), Biomass(2,:)/len, '-x', 'MarkerSize', 10, 'LineWidth', 3, 'Color', 'g')



toc

% Rainfall in time plot
R = Initial_rain;

J = ceil(abs(R-min_R)/Rainfall_drop);
 
% Increasing R from 1.25 to 1.4
J2 = ceil(abs(R-max_R)/Rainfall_drop);

time = zeros(1, J+J2+2);
rain = zeros(1, J+J2+2);

rain(1) = Initial_rain + Rainfall_drop*J2;


for i = 1:(J+J2)
    time(i+1) = time(i) + T;
    rain(i+1) = round(rain(i) - Rainfall_drop, 2);
end

time(end) = time(end-1)+T;
rain(end) = rain(end-1);

figure
stairs(time, rain, 'LineWidth', 3, 'Color', 'k')
set(gca, 'FontSize', 15)
xlabel('Time [day]', 'FontSize', 24)
ylabel('Rainfall [mm/day]',  'FontSize', 24)
ylim([0, max_R*1.01])
xlim([-100, time(end)*1.01])


% Sample plots
% Plot of Gap patterns at R = 1.20
figure
plot(xplot,PWO_IC_120(1,:), 'r', 'LineWidth', 3);
hold on
plot(xplot,PWO_IC_120(2,:), 'g', 'LineWidth', 3);
plot(xplot,PWO_IC_120(3,:), 'b', 'LineWidth', 3);
set(gca, 'FontSize', 15)
subtitle('R = 1.20',  'FontSize', 18);
xlabel("Spatial location [m]",  'FontSize', 18)
ylabel("Plant density [g/m^2], surface/soil water [mm]",  'FontSize', 18)
% legend('Plant biomass, P', 'Soil water, W', 'Surface water, O')
ylim([-0.1, 20])

% Plot of gap patterns at R = 0.8
figure
plot(xplot,PWO_IC_08(1,:), 'r', 'LineWidth', 3);
hold on
plot(xplot,PWO_IC_08(2,:), 'g', 'LineWidth', 3);
plot(xplot,PWO_IC_08(3,:), 'b', 'LineWidth', 3);
set(gca, 'FontSize', 15)
subtitle('R = 0.8',  'FontSize', 18);
xlabel("Spatial localtion [m]",  'FontSize', 18)
ylabel("Plant density [g/m^2], surface/soil water [mm]",  'FontSize', 18)
% legend('Plant biomass, P', 'Soil water, W', 'Surface water, O')
ylim([-0.1, 20])
