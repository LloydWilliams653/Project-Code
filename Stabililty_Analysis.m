%% Stability analysis
% clear
% Constants:
len = 100; % Domain size

syms P W O R k1 k2 c g d rW w0 alpha L

% Homogeneous equations:
dP_homog = c*g*W*P/(W+k1) - d*P;
dW_homog = alpha*O*(P+(k2*w0))/(P+k2) - g*W*P/(W+k1) - rW*W;
dO_homog = R - alpha*O*(P+(k2*w0))/(P+k2);

% Jacobian
a11P0 = simplify(diff(dP_homog, P,1)); % P = 0;
a11 = 0; % for P neq 0
a12 = simplify(diff(dP_homog, W,1));
a13 = simplify(diff(dP_homog, O,1));

a21 = simplify(diff(dW_homog, P,1));
a22 = simplify(diff(dW_homog, W,1));
a23 = simplify(diff(dW_homog, O,1));

a31 = simplify(diff(dO_homog, P, 1));
a32 = simplify(diff(dO_homog, W, 1));
a33 = simplify(diff(dO_homog, O, 1));

% Equilibria:
% For P neq 0:
W_Pneq = d*k1/(c*g - d);

P_neq = (R - rW*W)/(g*W) * (W+k1);
P_Pneq = simplify(subs(P_neq, W, W_Pneq));

O_neq = R/alpha * (P+k2)/(P+k2*w0);
O_Pneq = subs(O_neq, P, P_Pneq);

% For P = 0:
P_P0 = 0;
W_P0 = R/rW;
O_P0 = R/(alpha*w0);

% Jacobian Matrix for P neq 0
A = [a11, a12, a13 ;...
     a21, a22, a23; ...
     a31, a32, a33];

% Jacobian Matrix for P = 0
A_0 = [a11P0, 0, a13 ;...
     a21, a22, a23; ...
     a31, a32, a33];


% Matrices in terms of parameters & R only
A_Pneq_old = subs(A, [P, W, O], [P_Pneq, W_Pneq, O_Pneq]);

A_Pneq = simplify(subs(A_Pneq_old, [c, g, k1, k2, alpha, rW, w0, d], ...
    [10, 0.05, 5, 5, 0.2, 0.2, 0.2, 0.25]));

    

A_P0_old = simplify(subs(A_0, [P, W, O], ...
    [P_P0, W_P0, O_P0]));

A_P0 = simplify(subs(A_P0_old, [c, g, k1, k2, alpha, rW, w0, d], ...
    [10, 0.05, 5, 5, 0.2, 0.2, 0.2, 0.25]));


% Diffusion contants
DP = 0.1;
DW = 0.1;
DO = 100;

% Jacobian with diffusive terms

Jac_simp = [A_P0(1,1)-L^2*DP , 0 ; A_P0(2,1)+(A_P0(2,3)*A_P0(3,1))/((DO*L^2)-A_P0(3,3)), A_P0(2,2) - L^2*DW]; % For P = 0

Jac_simp_Pneq = [-L^2*DP , A_Pneq(1,2) ; A_Pneq(2,1)+(A_Pneq(2,3)*A_Pneq(3,1))/((DO*L^2)-A_Pneq(3,3)), A_Pneq(2,2) - L^2*DW]; % For P neq 0











% Dispersion curve for P neq 0, R >= 1

% Find the growth rates for P neq 0
h11 = @(x,y)max(real(eig(double(subs(Jac_simp_Pneq, [R,L], [x, y])))));

colors = {'r', 'g', 'b', 'm', 'c', 'k', 'y'};
i = 0; % for colors loop
num_functions = length(colors);
labels = cell(1, num_functions); 
lines = zeros(1, num_functions); 


figure
% Dispersion curves
for precip = 1:(1.3-1)/6:1.3
    i = i + 1;
    h = fplot(@(x)h11(precip,x), [0,1], 'LineWidth', 2, 'Color', colors{i}, 'DisplayName', sprintf(['R = ' num2str(precip)], i));
    labels{i} = sprintf(['R = ' num2str(precip)], i);
    lines(i) = h; 
    hold on
    yline(0)
    xline(0)
    set(gca,'FontSize', 12);
    xlabel('Spatial Frequency, L [m^{-1}]', 'FontSize', 24)
    ylabel('Growth rate, Re(\sigma)', 'FontSize', 24)
    % title('Dispersion Curves for P \neq 0', 'FontSize', 20)
end

legend(lines, labels, 'Location', 'best', 'Color', 'none');








% Minimimum eigenvalues:
% h2 = @(x,y)min(real(eig(double(subs(Jac_simp_Pneq, [R,L], [x, y]))))); %

% Max eigenvalues:
h1 = @(x, y) arrayfun(@(x_val, y_val) max(real(eig(double(subs(Jac_simp_Pneq, [R, L], [x_val, y_val]))))), x, y);

inc = 1200; % Number of spacings 

x = 1:(1.3-1)/inc:1.3; % Rainfall
y = 0:(0.7-0)/inc:0.7; % Spatial frequency
instab_ = zeros(length(x), length(y));
axes = zeros(2, length(y));


for i = 1:length(x)
    precip = x(i)*ones(1,length(x));
    axes(:,i) = [x(i); y(i)];
    instab_(i,:) = h1(precip, y);
end
instab = instab_ <= 0;


% Find the maximum rainfall for which patterns can form 
for i = 2:length(instab)
    if sum(instab(i-1,:)) ~= length(instab) && sum(instab(i,:)) == length(instab)
        max_unstab_R = axes(1,i); 
        disp(max_unstab_R)
        disp(i)
    end
end










%%
% Rainfall plot
figure
imagesc(axes(1,:), axes(2,:), instab')
axis xy
custom_map = [0, 0.5, 0; 1, 1, 1];  
colormap(custom_map);
set(gca, 'FontSize', 12)
ylabel('Spatial Frequency, L [m^{-1}]', 'FontSize', 24)
xlabel('Rainfall, R [mm/day]', 'FontSize', 24)
xlim([axes(1,1), axes(1,end)])
ylim([axes(2,1), axes(2,end)])


toc


j = 0;
colors = {'r', 'g', 'b', 'm', 'c', 'k', 'y'};
for i = axes(1,1):(axes(1,end)-axes(1,1))/6:axes(1,end)
    j = j + 1;
xline(i, 'Color', colors{j}, 'LineStyle', '--', 'LineWidth', 2)
end

%% Frequency lines plot
figure
imagesc(axes(1,:), axes(2,:), instab')
axis xy
custom_map = [0, 0.5, 0; 1, 1, 1];  % Dark green and white
colormap(custom_map);


spacings = 0:length(axes(2,1):2*pi/len:axes(2,end));
j = 0;
gradient_colormap = [linspace(0,1,length(spacings))', zeros(length(spacings),1), linspace(1,0,length(spacings))'];

for i = axes(2,1):2*pi/len:axes(2,end)
    j = j + 1;
yline(i, 'Color', gradient_colormap(j,:), 'LineWidth', 2, 'Label', sprintf('n = %d', j-1))
end
set(gca, 'FontSize', 12)

ylabel('Spatial Frequency, L [m^{-1}]', 'FontSize', 24)
xlabel('Rainfall, R [mm/day]', 'FontSize', 24)


