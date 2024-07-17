function [x] = Algo1(u, alpha, zeta, I, L, N, kron_inv_Z_alpha, A_p)
%ALGO1 Summary of this function goes here
%   Detailed explanation goes here
    
    phi =  2* pi * rand(L * N, 1);
    gamma = exp(1j *2* pi * rand(L * N, 1));
    x =  exp(1j *2* pi * rand(L * N, 1));

    i = 1;
    while i < I &&  norm(x - exp(1i * phi),2) > zeta
        b =  x + gamma; % b(i) = x(i) + gamma(i)
        phi = angle(b); % eq 20: phi(i+1) = angle (phi(i))
        x_hat = exp(1i * phi) - gamma; % x_hat (i) = exp(1i * phi(i+1)) - gamma(i);
        q = q_(u, x_hat, alpha, A_p); % q(i) = sum(Ap_hp,2) + x_hat(i); 
        x = kron_inv_Z_alpha * q; % eq22:  x(i+1) = kron(I_N,inv(Z_alpha)) * q(i);
        
        gamma = gamma + x - exp(1i * phi);  %18c: gamma (i+1) = gamma(i) + x(i+1) - exp(1i * phi(i+1);
        i = i+1;
        % norm(x - exp(1i * phi),2);
    end
end

