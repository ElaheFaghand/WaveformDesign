%% Stoica 2008 CA (cyclic algorithm) Method for Waveform design in article: 
%% Waveform Synthesis for Diversity-Based Transmit Beampattern Design
%% It is written by Elahe Faghand
clc
clear
tic
theta = (-90:.1:90)'*pi/180; % radian 181 points
N = 16; % ANTENNA NUMBER
L = 10; % SAMPLE NUMBER
%%
R = importdata("R.mat"); % Enter the saved R from SQP_matching or SQP_SLL or SQP_SLL_wideband 
P_R = Beam_Pattern (N, theta, R);
plot(theta*180/pi ,10*log10(abs(P_R))),hold on,  grid on
%% 
U = 0.7 * randn (L,N) + 0.7*1j * randn (L,N);
U_hat = ones (L,N);
c = sqrt(diag(R(1)));
for i = 1:1000
    % step 1
    X =  c * exp(1j * angle (sqrt(L) * U * sqrtm(R)));
    % step 2
    U_bar_sigma_U_tild_c = sqrt(L) * sqrtm(R) * X';
    [ U_bar , sigma , U_tild] = svd( U_bar_sigma_U_tild_c , 'econ');
    U_hat = U_tild * U_bar';
    i = i + 1;
    sum(sum (abs(U_hat - U)));
    U = U_hat;
end
% X = sqrt(L) * U_hat * sqrtm(R); % in this case beampattern is perfectly match but it is not CM
R_CA = 1/L * X' * X;
P_CA = Beam_Pattern (N, theta, R_CA);
plot(theta*180/pi ,10*log10(abs(P_CA))),hold on,  grid on
toc
