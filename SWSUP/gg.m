clc
clear
format long

K = 2001;
M = 20;
p_d = zeros(K, 1);

gamma_L = -0.5;
gamma_R = 0.5;

gamma = reshape(linspace(-1,1, K), K, 1);
p_d((gamma>gamma_L)&(gamma<gamma_R)) = 1;

% figure, hold on, plot(gamma, (p_d-1)*20)

eta_1 = gamma_L + 1/M;
eta_I = gamma_R - 1/M;

I = ceil((gamma_R-gamma_L)*M/2);

h = reshape(linspace(eta_1, eta_I, I), 1, I);
b = ones(1, I);

% disp(h)

%% update beta
for j =1:15
    x = gamma-h;
    a = pi/2;
    ax = a*x;
    Max = M*ax;
    
    P = 1/M.^2 * (sin(Max)./sin(ax)).^2;
    P(isnan(P)) = 1;
    
    
    p_l = P*b';
    g = -2/K * (p_d- p_l)' * P;
    H = 2/K * (P' * P);
    
    b = b - (g * inv(H));
    % b = b - (g / H);  % according to matlab, this is equivalent.
   
    % if j == 1
    %     plot(gamma, 10*log10(p_l))
    % end

    %% update eta
    for i=1:1000  
        x = gamma-h;
        a = pi/2;
        ax = a*x;
        Max = M*ax;
        
        sin_ax = sin(ax);
        sin_Max = sin(Max);
        sin2_ax = sin_ax.^2;
        sin2_Max = sin_Max .^ 2;
        
        cos_ax = cos(ax);
        cos_Max = cos(Max);
        cos2_ax = cos_ax.^2;
        cos2_Max = cos_Max.^2;
        
        P = 1/M^2 * (sin(Max)./sin_ax).^2;
        P(isnan(P)) = 1;

        Pdot = 2 * a * ((sin2_Max .* cos_ax) ./ sin_ax.^3 - M .* (sin_Max .* cos_Max) ./ sin2_ax);
        Pdot(isnan(Pdot)) = 1;

        Pdotdot = 2*a.^2 ./ sin_ax.^4 .* (M.^2 * (cos2_Max .* sin2_ax - sin2_Max .* sin2_ax) - 4 * M * sin_Max .* cos_Max .* sin_ax .* cos_ax + sin2_Max .* (sin2_ax + 3 * cos2_ax));
        Pdotdot(isnan(Pdotdot)) = 1;

        p_l = P * b';
        Edot = - 2/K * b .* ((p_d-p_l)' * Pdot);
        Edotdot = 2/K .* (b.^2 .* (ones(1, K) * Pdot.^2) - b.*((p_d-p_l)' * Pdotdot));
        
        h = h - (Edot./Edotdot);

    end

    MSE = 1/K * ((p_d - p_l)' * (p_d - p_l))
end

disp(h)
MSE = 1/K * ((p_d - p_l)' * (p_d - p_l));
disp(MSE)
plot(gamma, 10*log10(p_l))

