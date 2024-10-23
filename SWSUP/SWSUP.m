% Main Final version of SWSUP method in article:
% A Novel Method Based on Sequential Unconstrained Programming for Transmit Beamforming in Colocated MIMO Radars
% It is written by Elahe Faghand
clc
clear    
M = 30; % antenna number
gamma = (-1:0.001:1);  
theta = asin(gamma) * 180/pi;
    
%% Base covariance matrix construction
R_0 = ones(M); 
R_0 = R_0 / sum(sum(R_0));
       
%% eta values
% % single lobe
theta_start = -30; % degree
theta_stop = 30;  % degree
% 
gamma_L = sind(theta_start);
gamma_R = sind(theta_stop); 
P_d = (1.*(gamma > gamma_L) .* (gamma < gamma_R));

% P_d = [0.1*ones(1,200),zeros(1,200),0.1*ones(1,200),ones(1,600),0.1*ones(1,601)];
% % multi lobes
% theta_start = -15; % degree
% theta_stop = 15;  % degree
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

%% Eta definition that used in only beta variables
eta_1 = (gamma_L) + (1./M);
eta_I = (gamma_R) - (1./M);
I =  ceil((eta_I - eta_1) / (2/M)) + 1;
eta = linspace(eta_1, eta_I,I);

% % multi lobes
% eta_1 = (gamma_L2) + (1./M);
% eta_I = (gamma_R2) - (1./M);
% I =  ceil((eta_I - eta_1) / (2/M)) + 1;
% eta2 = linspace(eta_1, eta_I,I);

% eta_1 = (gamma_L3) + (1./M);
% eta_I = (gamma_R3) - (1./M);
% I =  ceil((eta_I - eta_1) / (2/M)) + 1;
% eta3 = linspace(eta_1, eta_I,I);

% eta = [eta eta2 eta3];

%%
beta_initial = 0.1 * ones(1, length(eta));
V = 10; % set V = 1 if you want only run weight evaluation
tic
for Iteration = 1:V
    %% Beta calculation
    [P_R_total1, beta, P_shift] = Beta_calculation (gamma, P_d, M, eta, beta_initial); % weight evaluation
    %% Theta calculation, set Iteration = 1:10 if you want  run beta eta code % shift evaluation
    [eta_updated, P_R_total] = Theta_calculation(gamma,  P_d, M, eta, beta);
    beta_initial = beta';
    eta = eta_updated;
end
toc
plot(theta , 10*log10(abs(P_R_total1))), hold on, grid on, axis([-90 90, -50 1])
    
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
