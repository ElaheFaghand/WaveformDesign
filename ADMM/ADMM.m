clc
clear all
tic
format long

gamma0 = (-1:0.01:1);  
theta = (asin(gamma0))' ; % radian

% theta = (-90:1:90)'*(pi/180); % radian
M = 16; % ANTENNA NUMBER
L = 128; % SAMPLE NUMBER
tic
%% one lobe
pl1 = sin(pi/180 * -17);
ph1 = sin(pi/180 * 17);
d = (1.*(sin(theta)>pl1).*(sin(theta)<ph1));

%% 3 lobes
% p=sin(pi/180*-50);
% ph=sin(pi/180*-40);
% p2=sin(pi/180*-15);
% ph2=sin(pi/180*15); 
% p3=sin(pi/180*40);
% ph3=sin(pi/180*60);
% d=((1.*(sin(theta)>p).*(sin(theta)<ph)))+(1.*(sin(theta)>p2).*(sin(theta)<ph2))+(1.*(sin(theta)>p3).*(sin(theta)<ph3));
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
while m < 60 % 36, 37, 38
    u = 1/rho_1 * u;
    v = 1/ rho_2 * v;
    
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
    gamma = rho_1/2 * (r - u ) + rho_2/2 * (T_prime)' * (1 - v ); %  (eq. 32)
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
    zeta = rho_1/2 * (h + u ) + rho_2/2 * (T)' * (1 - v ); %  (eq. 32)
    r = inv(omega) * zeta; % (eq. 33)

   %% update_u,v
    u = u + h - r;

    % update_v
    for i =2: (M*L+1)
        E = zeros (M*L + 1 , M*L + 1);
        E(i,i) = 1;
        T_prime_1 (:,i-1) = [E * r ]; % below (eq. 21)
    end
    T_prime = T_prime_1';

    v = v + T_prime * h - 1;


    %% MSE
%     for k = 1:length(theta)
%         obj (k) = abs(r'*A(:,:,k)*r)^2;
%     end
% obj = sum(obj)
%     % r = r/r(1);
% 
%     r1 = r/r(1);
%     r1(1)=[];
%     for k = 1:length(theta)
%         P (k) = (r1'*R(:,:,k)*r1);
%     end
% plot(theta*180/pi,(abs(P))),hold on,  grid on

%% MSE
   % h_r = norm(h - r);
   % G = norm(T*r-1);

        r1 = r;
        r1(1)=[];
        for k = 1:length(theta)
            P (k) = (r1'*R(:,:,k)*r1)/ abs(r(1))^2;
            err (k) = abs( P (k)-d(k))^2;
        end

        err = sum(err)/length(theta)
%   plot(theta*180/pi,(abs(P))),hold on,  grid on
  max_r(m) = max (abs(r1)); 
  min_r (m) = min (abs(r1));
  m = m + 1

end

% plot
plot(theta*180/pi,(10*log10(abs(P)))),hold on,  grid on
 figure
plot((1:m-1),  max_r), hold on,plot((1:m-1),  min_r)
MSE  = 1/(length(theta)) * ((d - P')' * (d - P'))
toc