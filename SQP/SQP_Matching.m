%% Stoica 2007 SQP Method for solving matching design problem in section III-C in article: 
%% On Probing Signal Design For MIMO Radar
%% It is written by Elahe Faghand
clc
clear all
format long
gamma = (-1:0.001:1);  % 2001 points are considered
% theta = (-90:0.1:90)' *pi/180;
theta = (asin(gamma))' ;
tt = theta * 180/pi;
L = length(theta);
t = sin(theta);
N = 30; % number of antennas
%% multi lobe desired beampattern definition
% p=sin(pi/180*-50);
% ph=sin(pi/180*-40);
% p2=sin(pi/180*-15);
% ph2=sin(pi/180*15);  
% p3=sin(pi/180*40);
% ph3=sin(pi/180*60);
% r=((1.*(t>p).*(t<ph)))+(1.*(t>p2).*(t<ph2))+(1.*(t>p3).*(t<ph3));
%% single lobe desired beampattern definition
theta_start = -30; % degree
theta_stop = 30;  % degree
pl1 = sin(pi/180 * theta_start);
ph1 = sin(pi/180 * theta_stop);
r=(1.*(t>pl1).*(t<ph1));
%% SQP Method
A = zeros(L,N^2);
for l = 1:L
    aT = exp(  1j * pi*(0:N-1)*sin(theta(l)) ).'; % steering vector
    A(l,:) = kron(aT.',aT');
end
cvx_begin quiet
variable R(N,N) complex  % covariance matrix
variable alphaa 
tic
minimize(norm(alphaa*r-1.*A*vec(R),2))
subject to
    R(logical(eye(N))) == 1/N;
    R == hermitian_semidefinite(N);
cvx_end
R = R/alphaa;
toc
for l = 1:L
    et = exp(1j*pi*(0:N-1)*sin(theta(l))).';
    P_sqp(l) = et'*(R)*et;
end
% plot(tt,(abs(P_sqp))),hold on,  grid on
 figure
 plot(tt,10*log10(abs(P_sqp)/max(abs(P_sqp))))
MSE = 1/L * ((r' - P_sqp) * (r' - P_sqp)')
