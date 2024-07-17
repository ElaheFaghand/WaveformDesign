% 9/1/2024 algorithm 1 started
% 12/5/2024 algorithm 1 is done :-)
% 11/06/2024 algorithm 2 started, I hope i can
clc
clear all
tic
rng(7)
format long
theta = (-90:1:90)'*(pi/180); % radian
L = 16; % ANTENNA NUMBER
N = 32; % SAMPLE NUMBER
zeta = 1e-8;
I = 500; % I = 2500;
alpha = 200;
rho = 0.01; % step size
T_0 = 5000;     
iota = 0; % i did not find the value in the text
t = 1;
%% desired beampattern
Theta = (-20:1:20)* (pi/180); % mainlobe region radian
Omega = [(-90:1:-28),(28:1:90)] * (pi/180); % sideloberegion radian
M = length(Theta);
S = length(Omega);
P = M + S;
%%
lambda_s = 1/sqrt(2) * (rand(N,S) + 1i* rand(N,S));
nu_m = 1/sqrt(2) * (rand(N,M) + 1i* rand(N,M));
I_N = eye(N);
I_L = eye(L);
% phi = rand(L*N,1);
% gamma = 1/sqrt(2) * (rand(L*N,1) + 1i* rand(L*N,1));
% x =  1/sqrt(2) * (rand(L*N,1) + 1i* rand(L*N,1));
for l = 1:length(sin(theta))
     aT = exp(  1j * pi*(0:L-1)*sin(theta(l)) ).';
     A_theta(:,:,l) = kron (I_N, aT);
end
eps_f  = [1 2];

while t <= T_0 %&&  abs(eps_f(t+1) - eps_f(t)) > iota
phi = rand(L*N,1);
gamma = 1/sqrt(2) * (rand(L*N,1) + 1i* rand(L*N,1));
x =  exp(1j *2*pi* phi);

%% Steering Vector
    for s = 1:S
        aTs(:,s) = exp( 1j * pi*(0:L-1)*sin(Omega(s)) ).'; % L * 1
        A_s (:,:,s) = kron(I_N,aTs(:,s)); % LN * N in each angle 
        y_s (:,s) = A_s(:,:,s)' * x; % N * 1
        u_s (:,s) = y_s (:,s) + lambda_s(:,s);
        aahs(:,:,s) = aTs(:,s) * aTs(:,s)';
        As_hs(:,s) = A_s(:,:,s) * u_s (:,s);
    end
    aah_s = sum (aahs,3);
    Ashs = sum(As_hs,2);

    for m = 1:M
        aTm(:,m) = exp( 1j * pi*(0:L-1)*sin(Theta(m)) ).';
        A_m (:,:,m) = kron(I_N,aTm(:,m));
        z_m (:,m) = A_m(:,:,m)' * x;
        u_m (:,m) = z_m (:,m) + nu_m(:,m);
        aahm(:,:,m) = aTm(:,m) * aTm(:,m)';
        Am_hm(:,m) =  A_m(:,:,m) * u_m (:,m);
    end

    aah_m = sum (aahm,3);
    Amhm = sum (Am_hm,2);

    sum_aah = aah_m + aah_s;
    sum_Ap_hp = Ashs + Amhm;

%% concatenate

%     u_p = [u_m, u_s]; % h_p = zeros(N,1);
%     aTp = [aTm, aTs]; 
%     A_p = cat(3,A_s,A_m);
%     for p = 1:P
%         aah(:,:,p) = aTp(:,p) * aTp(:,p)';
%         Ap_hp(:,p) =  A_p(:,:,p) * u_p (:,p);
%     end
  
%     sum_aah = sum (aah,3);
%     sum_Ap_hp = sum(Ap_hp,2);
%%
    Z_alpha = sum_aah + alpha/2 * I_L; 

%% Alg 1
i = 0;
    while i < I &&  norm(x - exp(1i * phi),2) > zeta
        b =  x + gamma; % b(i) = x(i) + gamma(i)
        phi = angle(b); % eq 20: phi(i+1) = angle (phi(i))
        x_hat = exp(1i * phi) - gamma; % x_hat (i) = exp(1i * phi(i+1)) - gamma(i);
        q = sum_Ap_hp + 0.5*alpha*x_hat; % q(i) = sum(Ap_hp,2) + x_hat(i); 
        x = kron(I_N , inv(Z_alpha)) * q; % eq22:  x(i+1) = kron(I_N,inv(Z_alpha)) * q(i);
        %   x = inv(B) * q; % eq22:  x(i+1) = kron(I_N,inv(Z_alpha)) * q(i);
       
        gamma = gamma + x - exp(1i * phi);  %18c: gamma (i+1) = gamma(i) + x(i+1) - exp(1i * phi(i+1);
        i = i+1;
        norm(x - exp(1i * phi),2)
    end
%% Alg 2 (13b)

    for s =1:S
        y_hat_s (:,s) = A_s(:,:,s)'*x - lambda_s(:,s); % below (23)
        norm_y_hat_s (s) =  norm(y_hat_s (:,s),2);
    end
    eta_hat = unique(round(norm_y_hat_s,2),'first');
    K = length (eta_hat)-1;
    for k = 2:K+1
        a = rho * (K- k +2);
        b = -rho * sum(eta_hat(k:K+1)) ;
        c = 2;
        delta = b^2 - 4*a*c;
        v1 = (-b - sqrt(delta))/(2*a);
        v2 = (-b + sqrt(delta))/(2*a);
    
        if  (delta <= 0) || (delta >= 0) %&& (v1 <= v2) && ( v2 <= eta_hat(k-1)) && (eta_hat(k-1) <=eta_hat(k)) )
            mu = eta_hat(k-1);
         disp ("case111111111")
        end
        if (delta >= 0) && (v1 <= eta_hat(k-1)) && ( eta_hat(k-1) <= eta_hat(k)) && (eta_hat(k) <=v2)
            mu = eta_hat(k);
            disp ("case2")
        end
        if  (delta >= 0) && (v1 <= eta_hat(k-1)) && ( eta_hat(k-1) <= v2) && (v2 <= eta_hat(k))
            mu = v2;
            disp ("case3333333333")
        end
            
        if (delta >= 0) && (eta_hat(k-1) <= v1) && (v1 <=v2) && (v2<= eta_hat(k))
            f1 = f_k(eta_hat(k-1), eta_hat, rho, k, K);
            f2 = f_k(v2, eta_hat, rho, k, K);
            if f1 < f2
                mu = eta_hat(k-1);
            else
                mu = v2;
            end
            disp ("case44444444444")
        end
        if (delta >= 0) && (eta_hat(k-1) <= v1) && (v1 <=eta_hat(k)) && (eta_hat(k) <= v2)
            
            f1 = f_k(eta_hat(k-1), eta_hat, rho, k, K);
            f3 = f_k(eta_hat(k), eta_hat, rho, k, K);
            
            if f1 < f3
                mu = eta_hat(k-1);
            else
                mu = eta_hat(k);
            end
            disp ("case5555555555")
        end
      mu_s(k-1) = mu;
      f(k-1) = f_k(mu, eta_hat, rho, k, K);
    end
    k_selct = find(f==min(f));
    eta = mu_s(k_selct(1))^2;
    
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
        norm_z_hat_m (m) =  norm(z_hat_m (:,m),2);
    end
    epsilon = unique(round(norm_z_hat_m,2),'first');
    K = length (epsilon)-1;
    for k = 2:K+1
        a_bar = rho * (k-1);
        b_bar = -rho * sum(epsilon(2:k)) ;
        c_bar = -2;
        g_p = g_prime(epsilon(k), a_bar, b_bar, c_bar); 
        
        delta_bar = b_bar^2 - 4*a_bar*c_bar;
        varsigma_1 = (-b_bar + sqrt(delta_bar))/(2*a_bar);
    
        a_bar_p = rho * (k-2);
        b_bar_p = -rho * sum(epsilon(2:k-1)) ;
        c_bar_p = -2;    
        g_p_p = g_prime(epsilon(k-1), a_bar_p, b_bar_p, c_bar_p); 
    
        if  (g_p_p > 0) 
            omega = epsilon(k-1);
        end
        if (g_p < 0) 
            omega = epsilon(k);
        end
        if  (g_p_p < 0) && (g_p > 0) 
            omega = varsigma_1;
        end
        
      omega_s(k-1) = omega;
      g(k-1) =  g_prime(omega, a_bar, b_bar, c_bar);
    end
    k_selct = find(g==min(g));
    eps_m = omega_s(k_selct(1))^2;
    
    for m =1:M
        if norm_z_hat_m(m) > sqrt(eps_m)
        z_m(:,m) = z_hat_m(:,m) ;
        else
        z_m(:,m) = sqrt(eps_m) * z_hat_m(:,m)/norm_z_hat_m(m);
       end
    end
    % 13d, 13e
    for s =1:S
        lambda_s(:,s) = lambda_s(:,s) + y_s(:,s) - A_s(:,:,s)' * x;
    end
    for m =1:M
        nu_m(:,m) = nu_m (:,m)+ z_m(:,m) - A_m(:,:,m)' * x;
    end
    t = t + 1;
    eps_f(t+1) = eps_m;
    plot(mu_s,f),grid on, hold on
    clear mu_s f omega_s g
   
-log10(eps_m/eta) ;
end
for l = 1:length(theta)
    P(l) = x'*A(:,:,l)*A(:,:,l)'*x;
end
plot(theta*180/pi,(10*log10(abs(P)))); grid on