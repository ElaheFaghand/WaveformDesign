% Final version of SWSUP method 
clc
clear 
i=1;
% for i =1: 1000
M = 30; % antenna number
gamma = (-1:0.001:1);  
theta = asin(gamma) * 180/pi;
%  theta = (-90:.1:90);
%  gamma = sind(theta);
%% Base covariance matrix construction
R_0 = ones(M); 
R_0 = R_0 / sum(sum(R_0));
%% eta values

% % single lobe
% a = rand;
% theta_start = randi([-75,-2]); % degree
% theta_stop = randi([2,75]);  % degree
% 
theta_start = -30; % degree
theta_stop = 30;  % degree
gamma_L = sind(theta_start);
gamma_R = sind(theta_stop); 

P_d = (1.*(gamma > gamma_L) .* (gamma < gamma_R));
% P_d = [0.1*ones(1,200),zeros(1,200),0.2*ones(1,200),ones(1,600),0.2*ones(1,601)];
% P_d = [zeros(1,600),ones(1,250),zeros(1,100),ones(1,250),zeros(1,601)];
% P_d = [zeros(1,600),ones(1,270),zeros(1,60),ones(1,270),zeros(1,601)];

% plot (theta,P_d)
% multi lobes
%%
% Multi lobe 
% theta_start = randi([0,40]);
% theta_stop = randi ([theta_start+5, 80]);
% gamma_L = sind(theta_start);
% gamma_R = sind(theta_stop); 
% theta_stop = randi([-40,0]);
% theta_start = randi ([-80,theta_stop-5]);
% gamma_L2 = sin(theta_start * pi/180);
% gamma_R2 = sin(theta_stop * pi/180);
% P_d = ((1.*(gamma>gamma_L).*(gamma<gamma_R)))+(1.*(gamma>gamma_L2).*(gamma<gamma_R2));
% figure
% plot (theta,P_d)
%%
% theta_start =5; % degree
% theta_stop = 5+l;  % degree
% 
% gamma_L2 = sin(theta_start * pi/180);
% gamma_R2 = sin(theta_stop * pi/180); 
% 
% theta_start = 40; % degree
% theta_stop = 60;  % degree
% 
% gamma_L3 = sin(theta_start * pi/180);
% gamma_R3 = sin(theta_stop * pi/180); 
% P_d = ((1.*(gamma>gamma_L).*(gamma<gamma_R)))+(1.*(gamma>gamma_L2).*(gamma<gamma_R2))+(1.*(gamma>gamma_L3).*(gamma<gamma_R3));
% P_d = ((1.*(gamma>gamma_L).*(gamma<gamma_R)))+(1.*(gamma>gamma_L2).*(gamma<gamma_R2));
% figure
% plot (theta,P_d)

%% Eta definition 1 for beta theta variables
% با این تعریف تغییر محسوس بین بتا خالی و زمانی که چندین ایتریشن از تتا
% میریم مشاهده میشه در حالی که با تعریف 3 هردو حالت یکی هستند ولی جوجه میگه
% چون MSE برامون مهمه بنابراین با تعریف 3 هست که به کمترین ارور ام اس ای
% میرسیم
% eta = [-2/M : -2/M : gamma_L, 2/M : 2/M:gamma_R, 0]
% I = length(eta);

%% Eta definition 3 that used in only beta variables
eta_1 = (gamma_L) + (1./M);
eta_I = (gamma_R) - (1./M);
I =  ceil((eta_I - eta_1) / (2/M)) + 1;
eta = linspace(eta_1, eta_I,I);

% % multi lobes
% eta_1 = (gamma_L2) + (1./M);
% eta_I = (gamma_R2) - (1./M);
% I =  ceil((eta_I - eta_1) / (2/M)) + 1;
% eta2 = linspace(eta_1, eta_I,I);
% 
% eta_1 = (gamma_L3) + (1./M);
% eta_I = (gamma_R3) - (1./M);
% I =  ceil((eta_I - eta_1) / (2/M)) + 1;
% eta3 = linspace(eta_1, eta_I,I);

% eta = [eta eta2];
% eta = [eta eta2 eta3];

%%
beta_initial = 0.1 * ones(1, length(eta));
V = 10; % set V = 1 if you want only run weight evaluation
tic
    for Iteration = 1:V
        %% Beta calculation
        [P_R_total1, beta, P_shift] = Beta_calculation (gamma, P_d, M, eta, beta_initial);
        %% Theta calculation, set Iteration = 1:10 if you want  run beta eta code
        [eta_updated, P_R_total] = Theta_calculation(gamma,  P_d, M, eta, beta);
        beta_initial = beta';
        eta = eta_updated;
%         MSE (Iteration) = (P_d - P_R_total')*(P_d' - P_R_total)/ length(theta);
    time(i) = toc;
    end
MSE (i) = (P_d - P_R_total')*(P_d' - P_R_total)/ length(theta);

figure
plot(theta , 10*log10(abs(P_R_total))), hold on, grid on, axis([-90 90, -50 1])
% hold on
% figure
% MSE(1) = .00394;
% plot(1:V,MSE);
% MSE(1) = .0129;
% MSE1 (i)= mean(MSE);
% i = i+1;
% clearvars -except i MSE time
% end
% mean(MSE)
% var (MSE)
% mean(time)
% var (time)
% plot (2*(4:4:60),MSE1)
%% Making R
% for i =1:length(eta)
%     % for i =1:length(eta_updated)
% 
%     R_L(:,:,i) = beta(i) * shift_R_construction (eta_updated(i),M,R_0);
% end
% R_L = sum(R_L,3);
% BP_L = Beam_Pattern (M,gamma,R_L);
% plot(gamma , 10*log10(abs(BP_L))), hold on, grid on %% plot only beta code
% R_L_shift = shift_R_construction (40 * pi/180,M,R_L);
% BP_L_shift = Beam_Pattern (M,gamma,R_L_shift);
% plot(gamma , 10*log10(abs(BP_L_shift))), hold on, grid on %% plot only beta code
% %%
% function R_shift = shift_R_construction (of,N,R)
%     A = zeros(N);
%         for i=1:N
%             A(i,i)= exp(1i*(pi)*(i-1)*(sin(of)));
%         end
%     R_shift=A*R*A';
% end
% 
% function Beam_Pattern = Beam_Pattern (N,gamma,R)
% Beam_Pattern = zeros(1, length(gamma));
%     for l = 1:length(gamma)
%         et = exp(1j*pi*(0:N-1)*(gamma(l))).';
%         Beam_Pattern(l) = et'*R*et;
%     end
% end