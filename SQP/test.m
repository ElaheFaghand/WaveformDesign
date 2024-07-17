clc
clear
N = 10;
M = 10;
theta = (-90:.1:90)'*pi/180; % radian 181 points
R_o = importdata("R.mat");
P_R = Beam_Pattern (N, theta, R_o);
plot(theta*180/pi ,10*log10(abs(P_R))),hold on,  grid on
delta = 1;
 while delta > 0.1   
    phi = rand(N,M);
    X = exp(1j*2*pi*phi);
    R = 1/N *(X'*X);
    delta = sum(sum(R - R_o))
 end
