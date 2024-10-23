% WBM method in article Constant Modulus MIMO Radar Waveform Design With Minimum Peak Sidelobe Transmit Beampattern
clc
clear all
    theta = (-90:1:90)'*(pi/180); % radian
    L = 16; % ANTENNA NUMBER
    N = 32; % SAMPLE NUMBER
    I_N = eye(N);
    zeta = 1e-8;
    I = 2500; % I = 2500;
    alpha = 200;
    rho = 0.005; % step size
    T_0 = 1000;  % 5000 is recommended in the article  
    iota = 0; % i did not find the value in the text
    t = 1;
    %% desired beampattern
    % Theta = [(-55:1:-35),(-10:1:10),(35:1:55)]* (pi/180); % mainlobe region radian
    % Omega = [(-90:1:-61),(-29:1:-16),(16:1:29),(61:1:90)] * (pi/180); % sideloberegion radian.
    
    Theta = (-20:1:20)* (pi/180); % mainlobe region radian
    Omega = [(-90:1:-28),(28:1:90)] * (pi/180); % sideloberegion radian
    % 
    ThetaOmega = [Theta,Omega];
    M = length(Theta);
    S = length(Omega);
    P = M + S;
    
    % x =  exp(1j *2* pi * rand(L * N, 1));
    lambda_s = rand(N,S);
    nu_m = rand(N,M);
    
    for s = 1:S
        a = exp( 1j * pi*(0:L-1)*sin(Omega(s)) ).'; % L * 1
        A_s (:,:,s) = kron(I_N, a); % LN * N in each angle 
        % y_s (:,s) = A_s(:,:,s)' * x; % N * 1
        % u_s (:,s) = y_s (:,s) + lambda_s(:,s);
    end
    for m = 1:M
        a = exp( 1j * pi*(0:L-1)*sin(Theta(m)) ).';
        A_m (:,:,m) = kron(I_N, a);
        % z_m (:,m) = A_m(:,:,m)' * x;
        % u_m (:,m) = z_m (:,m) + nu_m(:,m);
    end  
    I_L = eye(L);
    I_N = eye(N); 
    Z_alpha = (alpha / 2) * I_L;
    for  p = 1:length(ThetaOmega)
         a = exp( 1j * pi*(0 : L-1) * sin(ThetaOmega(p)) ).';
         A_p (:, :, p) = kron(I_N, a);
         Z_alpha = Z_alpha + a * a';
    end
    inv_Z_alpha = inv(Z_alpha);
    kron_inv_Z_alpha = kron(I_N , inv_Z_alpha);
    z_m  = rand(N, M);
    % nu_m = 1/sqrt(2) * (rand(N, M) + 1i * rand(N, M));
    
    y_s =  rand(N, S);
    % lambda_s = 1/sqrt(2) * (rand(N, S) + 1i* rand(N, S));
    
    eps_m = 1;
    eps_m_p = 0;
    tic
    while t <= T_0 &&  abs(eps_m - eps_m_p) > iota
        % delta = eps_m - eps_m_p
        eps_m_p = eps_m;
        u_m = z_m + nu_m;
        u_s = y_s + lambda_s;
        u = [u_m u_s];
        
        %% Alg 1
        x = Algo1(u, alpha, zeta, I, L, N, kron_inv_Z_alpha, A_p);
        
        %% Alg 2 (13b)
        for s =1:S
            y_hat_s (:,s) = A_s(:,:,s)'*x - lambda_s(:,s); % below (23)
            % norm_y_hat_s (s) =  norm(y_hat_s (:,s),2);
        end
        norm_y_hat_s  = vecnorm(y_hat_s,2);
        eta_hat = unique(round(vecnorm(y_hat_s,2),4),'first');
        K = length (eta_hat)-1;
        for k = 2:K+1
            a = rho * (K- k +2);
            b = -rho * sum(eta_hat(k:K+1)) ;
            c = 2;
            delta = b^2 - 4*a*c;
            v1 = (-b - sqrt(delta))/(2*a);
            v2 = (-b + sqrt(delta))/(2*a);
        
            if  (delta <= 0) ||  ((v1 <= v2) && ( v2 <= eta_hat(k-1)) && (eta_hat(k-1) <=eta_hat(k)) )|| ((eta_hat(k-1) <=eta_hat(k)) && (eta_hat(k) <= v1) && (v1<=v2)  )
                mu = eta_hat(k-1);
            end
            if (delta >= 0) && (v1 <= eta_hat(k-1)) && ( eta_hat(k-1) <= eta_hat(k)) && (eta_hat(k) <=v2)
                mu = eta_hat(k);
            end
            if  (delta >= 0) && (v1 <= eta_hat(k-1)) && ( eta_hat(k-1) <= v2) && (v2 <= eta_hat(k))
                mu = v2;
            end
                
            if (delta >= 0) && (eta_hat(k-1) <= v1) && (v1 <=v2) && (v2<= eta_hat(k))
                f1 = f_k(eta_hat(k-1), eta_hat, rho, k, K);
                f2 = f_k(v2, eta_hat, rho, k, K);
                if f1 < f2
                    mu = eta_hat(k-1);
                else
                    mu = v2;
                end
            end
            if (delta >= 0) && (eta_hat(k-1) <= v1) && (v1 <=eta_hat(k)) && (eta_hat(k) <= v2)
                f1 = f_k(eta_hat(k-1), eta_hat, rho, k, K);
                f3 = f_k(eta_hat(k), eta_hat, rho, k, K);
                
                if f1 < f3
                    mu = eta_hat(k-1);
                else
                    mu = eta_hat(k);
                end
            end
            mu_s(k-1) = mu;
            f(k-1) = f_k(mu, eta_hat, rho, k, K);
        end
        k_selct = find(f==min(f));
        eta = mu_s(k_selct(1))^2;
        %%
        for s =1:S
            if norm_y_hat_s(s) > sqrt(eta)
            y_s(:,s) = sqrt(eta) * y_hat_s(:,s)/norm_y_hat_s(s);
            else
            y_s(:,s) =  y_hat_s(:,s);
           end
        end
        
        % (13b)
        for m =1:M
            z_hat_m (:,m) = A_m(:,:,m)'*x - nu_m(:,m); % below (34)
        end
        norm_z_hat_m = vecnorm(z_hat_m, 2);
        epsilon = unique(round(norm_z_hat_m, 4),'first');
        K = length(epsilon)-1;
    
        % omega = zeros(K+1,1);
        % g_omega = zeros(K+1,1);
    
        for k = 2:K+1
            a_bar = rho * (k-1);
            b_bar = -rho * sum(epsilon(1:k-1)) ;
            c_bar = -2;
            % g_p = g_prime(epsilon(k), a_bar, b_bar, c_bar); 
            
            delta_bar = b_bar^2 - 4*a_bar*c_bar;
            varsigma_1 = (-b_bar + sqrt(delta_bar))/(2*a_bar);
                 
            if (g_prime(epsilon(k-1), a_bar, b_bar, c_bar) > 0) 
                omega(k) = epsilon(k-1);
            
            elseif (g_prime(epsilon(k), a_bar, b_bar, c_bar) < 0) 
                omega(k) = epsilon(k);
            
            elseif  ((g_prime(epsilon(k-1), a_bar, b_bar, c_bar) < 0) && (g_prime(epsilon(k), a_bar, b_bar, c_bar) > 0) )
                omega(k) = varsigma_1;
            else
                disp('ERROR')
            end
            
            g_omega(k) = g_k(omega(k), epsilon, rho, k);
        end
        omega_s = omega(2:end);
        g = g_omega(2:end);
        k_selct = find(g==min(g));
        eps_m = omega_s(k_selct)^2;
        
        for m =1:M
            if norm_z_hat_m(m) > sqrt(eps_m)
            z_m(:,m) = z_hat_m(:,m) ;
            else
            z_m(:,m) = sqrt(eps_m) * z_hat_m(:,m)./norm_z_hat_m(m);
           end
        end
        % 13d, 13e
        for s =1:S
            lambda_s(:,s) = lambda_s(:,s) + y_s(:,s) - A_s(:,:,s)' * x;
        end
        for m =1:M
            nu_m(:,m) = nu_m (:,m)+ z_m(:,m) - A_m(:,:,m)' * x;
        end
        t = t + 1
        % plot(mu_s,f),grid on, hold on
        % plot(omega_s,g),grid on, hold on
      % eps_m = eps_m
      % eta = eta
      % eps_m2eta = (eps_m/eta) 
      % -10*log10(eps_m/eta)
    end
toc
for l = 1:length(sin(theta))
     aT = exp(  1j * pi*(0:L-1)*sin(theta(l)) ).';
     A_theta(:,:,l) = kron (I_N, aT);
end
for l = 1:length(theta)
    P(l) = x'*A_theta(:,:,l)*A_theta(:,:,l)'*x;
end
% plot(theta*180/pi,(10*log10(abs(P)/max(abs(P))))); grid on
plot(theta*180/pi,(10*log10(abs(P))-37)); grid on
