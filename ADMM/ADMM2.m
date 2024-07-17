% It is proposed by Dr.merhshahi as considering rho1 and 
% rho2 as penalty functions and delet lagrangian (u,v) variables.

clc
clear all
tic
theta = (-90:1:90)'*(pi/180); % radian
M = 10; % ANTENNA NUMBER
L = 32; % SAMPLE NUMBER
tic
%% one lobe
pl1 = sin(pi/180*-30);
ph1 = sin(pi/180*30);
d = (1.*(sin(theta)>pl1).*(sin(theta)<ph1));

%% 3 lobes
% p=sin(pi/180*-50);
% ph=sin(pi/180*-30);
% p2=sin(pi/180*-10);
% ph2=sin(pi/180*10); 
% p3=sin(pi/180*30);
% ph3=sin(pi/180*50);
% d=((1.*(sin(theta)>p).*(sin(theta)<ph)))+(1.*(sin(theta)>p2).*(sin(theta)<ph2))+(1.*(sin(theta)>p3).*(sin(theta)<ph3));
% plot (theta*180/pi, d), hold on,

%% initial values
h = zeros(L*M +1 ,1) ; % h \in C (ML+1)
r = 1/sqrt(2) * (rand(L*M +1 ,1) + 1i* rand(L*M + 1,1)); % 19, h-r=0
% u = zeros(L*M +1 ,1); % u ~ y \in C (ML+1);
% v = zeros(L*M ,1); % v ~ z \in C (ML);
rho_1 = 2.1;
rho_2 = 2.1;
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
while m < 50 % 36, 37, 38
    %% update_h
        for i =2: (M*L+1)
            E = zeros (M*L + 1 , M*L + 1);
            E(i,i) = 1;
            T_prime_1 (:,i-1) = [E * r ]; % below (eq. 21)
        end
        T_prime = T_prime_1';
        
        for k = 1:length(theta)
             khi1 (:,:,k) = A (:,:,k) * r * r' * A (:,:,k)' ;
        end
    
        khi = (1/(length(theta))) * sum (khi1,3) + rho_1/2 * I  + rho_2/2 * T_prime' * T_prime; %  (eq. 31)
        gamma = rho_1/2 * (r  ) + rho_2/2 * (T_prime)' * (ones(320,1)); %  (eq. 32)
        h = inv(khi) * gamma; % (eq. 30)
    
        %% update_r
        for i =2: (M*L+1)
            E = zeros (M*L + 1 , M*L + 1);
            E(i,i) = 1;
            T (i-1,:) = [h' * E ]; % below (eq. 21)
        end
    
        for k = 1:length(theta)
             omega_1 (:,:,k) = A (:,:,k)' * h * h' * A (:,:,k) ;
        end
    
        omega = 1/(length(theta)) * sum (omega_1,3) + rho_1/2 * I + rho_2/2 * T' * T ; %  (eq. 31)
        zeta = rho_1/2 * (h ) + rho_2/2 * (T)' * ((ones(320,1))); %  (eq. 32)
        r = inv(omega) * zeta; % (eq. 33)
    
       %% update_u,v
    %     u = u + h - r;
    
        % update_v
    %     for i =2: (M*L+1)
    %         E = zeros (M*L + 1 , M*L + 1);
    %         E(i,i) = 1;
    %         T_prime_1 (:,i-1) = [E * r ]; % below (eq. 21)
    %     end
    %     T_prime = T_prime_1';
    % 
    %     v = v + T_prime * h - 1;
    
    rho_1 = 1.1*rho_1;
    rho_2 = 1.1*rho_2;
    
        %% MSE
   h_r = norm(h - r)
   G = norm(T*r-1)

        r1 = r;
        r1(1)=[];
        for k = 1:length(theta)
            P (k) = (r1'*R(:,:,k)*r1)/ abs(r(1))^2;
            err (k) = abs( P (k)-d(k))^2;
        end

        err = sum(err)/length(theta)

%         plot(theta*180/pi,(abs(P))),hold on,  grid on
  max_r(m) = max (abs(r1)); 
  min_r (m) = min (abs(r1));

  m = m + 1

end

% plot
plot(theta*180/pi,(abs(P))),hold on,  grid on
figure
plot((1:m-1),  max_r), hold on,plot((1:m-1),  min_r)
toc




