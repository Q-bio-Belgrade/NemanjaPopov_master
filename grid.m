%% 
clc;
close all;
clear all;
load perturbed_data_high.mat
plot(tExp,crRNA_exp,'bo');
hold on
%%
a = linspace(1,3,10);
b = linspace(0.2728,0.3222,10);
c = linspace(11.3367,16,10);
numCombinations = length(a) * length(b) * length(c);
array = zeros(3, numCombinations);
index = 1;

for i = 1:length(a)
    for j = 1:length(b)
        for k = 1:length(c)
            array(:, index) = [a(i); b(j); c(k)];
            index = index + 1;
        end
    end
end

sse_array = zeros(1, numCombinations);
opt_params_array = zeros(3,numCombinations);

parfor i = 1:numCombinations
        H0 = array(1, i);
        delta_H = array(2, i);
        n = array(3, i);
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
        tspan = 1:150;
        sse = SSE_FUNCTION([H0,delta_H,n], tExp,tspan, crRNA_exp, m, phi, k_star, lambda_pre, k, lambda, d, alpha, gamma, psi, K);
        
        sse_array(i) = sse;
        opt_params_array(:,i) = [H0,delta_H,n];
end
        
[min_sse, min_idx] = min(sse_array);
best_params = opt_params_array(:, min_idx);
disp(['Minimum SSE: ', num2str(min_sse)]);
disp(['Najbolji parametri: H0=', num2str(best_params(1)), ', delta_H=', num2str(best_params(2)), ', n=', num2str(best_params(3))]);

H0 = best_params(1);
delta_H = best_params(2);
n = best_params(3);
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

tspan = 0:149;
x0 = [Z_eq(end,1),Z_eq(end,2),Z_eq(end,3),Z_eq(end,4),Z_eq(end,5),Z_eq(end,6),Z_eq(end,7),Z_eq(end,8)];
[T_int, Z_int] = ode45(@(t, z) CRISPR(t, z, m, phi, k_star, lambda_pre, k, lambda, H0, delta_H, n, d, alpha, gamma, psi, K), tspan, x0);
crRNA_int = Z_int(:, 8);
plot(T_int,crRNA_int,'b-','LineWidth',2);
%%
function sse = SSE_FUNCTION(params, tExp,tspan, crRNA_exp, m, phi, k_star, lambda_pre, k, lambda, d, alpha, gamma, psi, K)
H0 = params(1);
delta_H = params(2);
n = params(3);

tspan_eq = 1:1500;
z0 = [10,10,10,10,10,10,10,10];
[T_eq, Z_eq] = ode23s(@(t, z_eq) CRISPR_eq(t, z_eq, m, phi, k_star, lambda_pre, k , lambda , H0 , n, d, alpha, gamma, psi, K), tspan_eq, z0);

x0 = [Z_eq(end,1),Z_eq(end,2),Z_eq(end,3),Z_eq(end,4),Z_eq(end,5),Z_eq(end,6),Z_eq(end,7),Z_eq(end,8)];
[T_int, Z_int] = ode45(@(t, z) CRISPR(t, z, m, phi, k_star, lambda_pre, k, lambda, H0, delta_H, n, d, alpha, gamma, psi, K), tspan, x0);
crRNA_int = Z_int(:, 8);

crRNA_fit = interp1(T_int,crRNA_int,tExp, 'linear', 'extrap');
sse = sum((crRNA_fit - crRNA_exp).^2);
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

%%
function zdot = CRISPR(t, z, m, phi, k_star, lambda_pre, k, lambda, H0,delta_H, n, d, alpha, gamma, psi, K)
H_hat = H0*(1-delta_H);
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