% It is the same as ADMM file, I only changed some parameters to reduce the
% time consumption for presenting the time in the articl while the
% constraints are the as close as to parametes article Jianli
clc
clear all

for ii = 1 : 10
    tic
    gamma0 = (-1:0.001:1);  
    theta = (asin(gamma0))' ; % radian
    K = length(theta);
    % theta = (-90:1:90)'*(pi/180); % radian
    M = 30; % ANTENNA NUMBER
    L = 10; % SAMPLE NUMBER
    E = zeros (M*L + 1 , M*L + 1);

    %% one lobe
% theta_start = randi([-75,-2]); % degree
% theta_stop = randi([2,75]);  % degree
%     pl1 = sin(pi/180 * theta_start);
%     ph1 = sin(pi/180 * theta_stop);
%     d = (1.*(sin(theta)>pl1).*(sin(theta)<ph1));
    %%
    theta_start1 = randi([0,40]);
    theta_stop1 = randi ([theta_start1+5, 80]);
    theta_stop = randi([-40,0]);
    theta_start = randi ([-80,theta_stop-5]);

p=sin(pi/180*theta_start1);
ph=sin(pi/180*theta_stop1);
p2=sin(pi/180*theta_start);
ph2=sin(pi/180*theta_stop);  
d=((1.*(sin(theta)>p).*(sin(theta)<ph)))+(1.*(sin(theta)>p2).*(sin(theta)<ph2));
       %% 3 lobes
%     p=sin(pi/180*-50);
%     ph=sin(pi/180*-40);
%     p2=sin(pi/180*-15);
%     ph2=sin(pi/180*15); 
%     p3=sin(pi/180*40);
%     ph3=sin(pi/180*60);
%     d = ((1.*(sin(theta)>p).*(sin(theta)<ph)))+(1.*(sin(theta)>p2).*(sin(theta)<ph2))+(1.*(sin(theta)>p3).*(sin(theta)<ph3));
    % plot (theta*180/pi, d)
    
    %% initial values
    h = zeros(L*M +1 ,1) ; % h \in C (ML+1)
    r = 1/sqrt(2) * (rand(L*M +1 ,1) + 1i* rand(L*M + 1,1)); % 19, h-r=0
    u = zeros(L*M +1 ,1); % u ~ y \in C (ML+1);
    v = zeros(L*M ,1); % v ~ z \in C (ML);
    rho_1 = 20;
    rho_2 = 20;
    eps_abs = 0.001;
    eps_rel = 0.01;
    
    I_L = eye(L);
    I = eye (M*L + 1);
    for k = 1:length(theta)
        aT = exp( - 1j * pi*(0 : M - 1)*sin(theta(k)) ).';
        R1 = kron(I_L,aT.');
        R(:,:,k) = R1'*R1;
        A (:,:,k) = [d(k),(zeros(1,M*L));(zeros(M*L,1)),-R(:,:,k)];
    end
    
    %%
       
    m = 1;
    while m < 15 % 36, 37, 38
        u = 1/rho_1 * u;
        v = 1/ rho_2 * v;
        
        %% update_h
        for i =2: (M*L+1)
            E(i,i) = 1;
            T_prime_1 (:,i-1) = [E * r ]; % below (eq. 21)
            E(i,i) = 0;
        end
        T_prime = T_prime_1';
        
        for k = 1:length(theta)
             khi1 (:,:,k) = A (:,:,k) * r * r' * A (:,:,k)' ;
        end
    
        khi = (1/K) * sum (khi1,3) + rho_1/2 * I  + rho_2/2 * T_prime' * T_prime; %  (eq. 31)
        gamma = rho_1/2 * (r - u ) + rho_2/2 * (T_prime)' * (1 - v ); %  (eq. 32)
        h = inv(khi) * gamma; % (eq. 30)
    
        %% update_r
        for i =2: (M*L+1)
            % E = zeros (M*L + 1 , M*L + 1);
            E(i,i) = 1;
            T (i-1,:) = [h' * E ]; % below (eq. 21)
            E(i,i) = 0;
        end
    
        for k = 1:length(theta)
             omega_1 (:,:,k) = A (:,:,k)' * h * h' * A (:,:,k) ;
        end
    
        omega = 1/K * sum (omega_1,3) + rho_1/2 * I + rho_2/2 * T' * T ; %  (eq. 31)
        zeta = rho_1/2 * (h + u ) + rho_2/2 * T' * (1 - v ); %  (eq. 32)
        r = inv(omega) * zeta; % (eq. 33)
    
       %% update_u,v
        u = u + h - r;
    
        % update_v
        for i =2: (M*L+1)
            % E = zeros (M*L + 1 , M*L + 1);
            E(i,i) = 1;
            T_prime_1 (:,i-1) = [E * r ]; % below (eq. 21)
            E(i,i) = 0;
    
        end
        T_prime = T_prime_1';
    
        v = v + T_prime * h - 1;  
    %% MSE
       % h_r = norm(h - r);
       % G = norm(T*r-1);
                r1 = r;
            r1(1)=[];
            for k = 1:length(theta)
                P (k) = (r1'*R(:,:,k)*r1)/ abs(r(1))^2;
            end
      % MSE  = 1/(length(theta)) * ((d - P')' * (d - P'))
               
      % plot(theta*180/pi,(abs(P))),hold on,  grid on
      max_r(m) = max (abs(r)); 
      min_r (m) = min (abs(r));

      m = m + 1;
   end
       r1 = r;
            r1(1)=[];
            for k = 1:length(theta)
                P (k) = (r1'*R(:,:,k)*r1)/ abs(r(1))^2;
            end
%             figure
%      plot(theta*180/pi,(10*log10(abs(P)))),hold on,  grid on
%     figure
%     plot((1:m-1),  max_r), hold on,plot((1:m-1),  min_r)
    MSE(ii) = 1/(length(theta)) * ((d - P')' * (d - P'));
    time (ii) = toc
    clearvars -except  ii MSE time
end
time = time/2.22;
mean(MSE)
var (MSE)
mean(time)
var (time)
toc