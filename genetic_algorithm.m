clc;
close all;
clear all;
load perturbed_data_high.mat
figure;
hold on

%%
lb = [1, 0.1, 5];
ub = [10, 0.3, 25];
matrix = zeros(200,3);

for i = 1:200
sigmaNoise = 15;
noise = normrnd(crRNA_exp,sigmaNoise);
x_pred = @(x) sum((predicted_x(tExp, x) - noise).^2);
options = optimoptions('ga','UseParallel',true);
optimalParams = ga(x_pred, 3, [], [], [], [], lb, ub,[],options);
matrix(i,:) = optimalParams;
opt_H0 = optimalParams(1);
opt_delta_H = optimalParams(2);
opt_n = optimalParams(3);
H0 = opt_H0;
delta_H = opt_delta_H;
n = opt_n;
phi = 10;
m = 2/3 * n;
k_star = 0.02;
lambda_pre = 1;
k = 3;
lambda = 1/5;
d = 0.0438;
alpha = 97.6174;
gamma = 55.9837;
psi = 57.6095;
K = 71.4272;
tspan_eq = 0:1500;
z0 = [10,10,10,10,10,10,10,10];
[T_eq, Z_eq] = ode23s(@(t, z_eq) CRISPR_eq(t, z_eq, m, phi, k_star, lambda_pre, k , lambda , H0 , n, d, alpha, gamma, psi, K), tspan_eq, z0);
tspan = 1:150;
x0 = [Z_eq(end,1),Z_eq(end,2),Z_eq(end,3),Z_eq(end,4),Z_eq(end,5),Z_eq(end,6),Z_eq(end,7),Z_eq(end,8)];     
[T_opt, Z_opt] = ode45(@(t, z) CRISPR(t, z, m, phi, k_star, lambda_pre, k, lambda, H0, delta_H, n, d, alpha, gamma, psi, K), tspan, x0);
crRNA_opt = Z_opt(:,8);
p1 = plot(T_opt,crRNA_opt,'b-',"LineWidth",1);
end

n = 18.3092;
H0 = 1.665;
delta_H = 0.2084;
phi = 10;
m = 2/3 * n;
k_star = 0.02;
lambda_pre = 1;
k = 3;
lambda = 1/5;
d = 0.0438;
alpha = 97.6174;
gamma = 55.9837;
psi = 57.6095;
K = 71.4272;

tspan_eq = 1:1500;
z0 = [10,10,10,10,10,10,10,10];
[T_eq, Z_eq] = ode23s(@(t, z_eq) CRISPR_eq( t , z_eq, m, phi, k_star, lambda_pre, k , lambda , H0 , n, d, alpha, gamma, psi, K), tspan_eq, z0);

tspan = 1:150;
x0 = [Z_eq(end,1),Z_eq(end,2),Z_eq(end,3),Z_eq(end,4),Z_eq(end,5),Z_eq(end,6),Z_eq(end,7),Z_eq(end,8)];
[T_sim, Z_sim] = ode45(@(t, z) CRISPR(t, z, m, phi, k_star, lambda_pre, k, lambda, H0, delta_H, n, d, alpha, gamma, psi, K), tspan, x0);
cr_sim = Z_sim(:, 8);
p2 = plot(tExp,crRNA_exp,'go','MarkerSize',5,'MarkerFaceColor','g');
p3 = plot(T_sim, cr_sim, 'r-','LineWidth',2);
set(gca,'Layer','top');
title('Genetski algoritam');
xlabel('Vreme(minuti)');
ylabel('Broj crRNK molekula');
legend([p1, p2, p3], {'Dinamika dobijena optimizacijom', 'Simulirani podaci', 'Originalna dinamika'}, 'Location', 'southeast');
hold off

figure(2);
[~,ax] = plotmatrix(matrix);
ylabel(ax(1,1),'H0');
ylabel(ax(2,1),'\DeltaH0');
ylabel(ax(3,1),'n');
xlabel(ax(3,1),'H0');
xlabel(ax(3,2),'\DeltaH0');
xlabel(ax(3,3),'n');
title('Matrica distribucije parametara');

figure(3);
tiledlayout(1,3);
nexttile;
boxplot(matrix(:,1));
yline(1.665,'--','Original');
ylabel('Ho');
nexttile;
boxplot(matrix(:,2));
yline(0.2084,'--','Original');
ylabel('\DeltaH0');
nexttile;
boxplot(matrix(:,3));
yline(18.3092,'--','Original');
ylabel('n');

%%
function zdot = CRISPR(t, z, m, phi, k_star, lambda_pre, k, lambda, H0, delta_H, n, d, alpha, gamma, psi, K)
    H_hat = H0 * (1 - delta_H);
    zdot = [
        alpha * ((1 + 2 * H_hat^m * (z(4) / K)) / (1 + H_hat^n + H_hat^m * (z(4) / K) + z(2) / (4 * 50 * K))) - lambda * z(1); 
        k * z(1) - d * z(2); 
        gamma * ((1 + z(2) / (4 * K)) / (1 + H_hat^n + z(2) / (4 * K))) - lambda * z(3); 
        k * z(3) - d * z(4); 
        psi * ((1 + z(2) / (4 * K)) / (1 + H_hat^n + z(2) / (4 * K))) - lambda * z(5); 
        k * z(5) - d * z(6); 
        phi - (lambda_pre + k_star * z(6)) * z(7); 
        k_star * z(6) * z(7) - d * z(8); 
    ];
end
%%
function x_pred = predicted_x(tExp, x)
        H0 = x(1);
        delta_H = x(2);
        n = x(3);
        phi = 10;
        m = 2/3 * n;
        k_star = 0.02;
        lambda_pre = 1;
        k = 3;
        lambda = 1/5;
        d = 0.0438;
        alpha = 97.6174;
        gamma = 55.9837;
        psi = 57.6095;
        K = 71.4272;

        tspan_eq = 1:1500;
        z0 = [10,10,10,10,10,10,10,10];
        [T_eq, Z_eq] = ode23s(@(t, z_eq) CRISPR_eq(t, z_eq, m, phi, k_star, lambda_pre, k , lambda , H0 , n, d, alpha, gamma, psi, K), tspan_eq, z0);

        tspan = 1:150;
        x0 = [Z_eq(end,1),Z_eq(end,2),Z_eq(end,3),Z_eq(end,4),Z_eq(end,5),Z_eq(end,6),Z_eq(end,7),Z_eq(end,8)];     
        [T_sim, Z_sim] = ode45(@(t, z) CRISPR(t, z, m, phi, k_star, lambda_pre, k, lambda, H0, delta_H, n, d, alpha, gamma, psi, K), tspan, x0);
        cr_sim = Z_sim(:, 8);
        x_pred = interp1(T_sim, cr_sim, tExp, 'linear', 'extrap');
end
%%
function z_eq = CRISPR_eq( ~ , z_eq, m, phi, k_star, lambda_pre, k , lambda , H0 , n, d, alpha, gamma, psi, K )
H_hat_eq = H0;
    z_eq = [
            alpha*((1+2*H_hat_eq^m*(z_eq(4)/K))/(1+H_hat_eq^n+H_hat_eq^m*(z_eq(4)/K)+z_eq(2)/(4*K*50))) - lambda*z_eq(1) ; 
            k*z_eq(1) - d*z_eq(2) ; 
            gamma*((1+z_eq(2)/(4*K))/(1+H_hat_eq^n+z_eq(2)/(4*K))) - lambda*z_eq(3) ; 
            k*z_eq(3) - d*z_eq(4) ; 
            psi*((1+z_eq(2)/(4*K))/(1+H_hat_eq^n+z_eq(2)/(4*K))) - lambda*z_eq(5) ; 
            k*z_eq(5) - d*z_eq(6); 
            phi - (lambda_pre + k_star*z_eq(6))*z_eq(7) ; 
            k_star*z_eq(6)*z_eq(7) - d*z_eq(8) ; 
        ];
end