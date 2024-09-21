clc;
close all;
clear all;
load perturbed_data_high.mat
plot(tExp,crRNA_exp,'go','MarkerSize',5,'MarkerFaceColor','g');
hold on
%%
phi = 10;
n = 18.3092;
m = 2/3 * n;
k_star = 0.02;
lambda_pre = 1;
k = 3;
lambda = 1/5;
H0 = 1.665;
delta_H = 0.2084;
d = 0.0438;
alpha = 97.6174;
gamma = 55.9837;
psi = 57.6095;
K = 71.4272;

z0 = [10,10,10,10,10,10,10,10];
tspan_eq = 0:1500;
[T_eq,Z_eq] = ode23s(@(t,z_eq)CRISPR_eq(t, z_eq, m, phi, k_star, lambda_pre, k , lambda , H0 , n, d, alpha, gamma, psi, K), tspan_eq, z0);


tspan = 0:149;
x0 = [Z_eq(end,1),Z_eq(end,2),Z_eq(end,3),Z_eq(end,4),Z_eq(end,5),Z_eq(end,6),Z_eq(end,7),Z_eq(end,8)];


[T, Z] = ode45(@(t, z) CRISPR(t, z, m, phi, k_star, lambda_pre, k, lambda, H0, delta_H, n, d, alpha, gamma, psi, K), tspan, x0);
crRNA = Z(:, 8);
plot(T, crRNA, 'r-','LineWidth',2);
xlabel('Vreme(minuti');
ylabel('Broj crRNK molekula');
legend('Simulirani podaci', 'Originalna dinamika','Location','southeast');
title('Dinamika crRNK molekula');

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