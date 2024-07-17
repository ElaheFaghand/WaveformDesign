clc
clear
theta = (-90:.1:90)'*pi/180; % radian 181 points
M = 16;
N = 1000000;
% gamma = 0.5;
% r = (0:M-1);
% R = (gamma.^abs(r-r'))
% P = Beam_Pattern (M, theta, R);
% plot(theta*180/pi ,10*log10(abs(P))),hold on,  grid on
%%
R = importdata("R.mat");
P = Beam_Pattern (M, theta, R);
plot(theta*180/pi ,10*log10(abs(P))),hold on,  grid on
[W,Lambda] = eig(R);
X = randn(N,M) * sqrtm(Lambda) * W';
R_N = X' * X/ N;
P_N = Beam_Pattern (M, theta, R_N);
plot(theta*180/pi ,10*log10(abs(P_N))),hold on,  grid on
figure
imagesc(abs(R))
figure
imagesc(abs(R_N))