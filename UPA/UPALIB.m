%% UPALIB method in article: 
%% Complexity Reduction in Beamforming of Uniform Array Antennas for MIMO Radars Section III.A
%% It is written by Elahe Faghand
clc
clear
format long
N = 30; % number of antennas
of = 1; 
i = 1;
M1 = -30; % strart point
M2 = 30;  % stop point
%%
gamma = (-1:0.001:1);  % 2001 points are considered
theta = (asin(gamma))' ;
tt = theta * 180/pi;
L = length(theta);
%%
m1 = -M2;
m2 = -M1;
%% making members of UPA library
for p = m1 : of : m2
o = exp(1i * pi * (0 : N-1) * tan(p * pi/180));
pp0(i,:,:) = o'*o;
i = i + 1;
end
m = m2 - m1 + 1;
t = sin(theta);
o = zeros(N);
%%  one lobe desired beam pattern
% pl1 = sin(pi/180 * M1);
% ph1 = sin(pi/180 * M2);
% r = (1.*(t>pl1).*(t<ph1));
%% 2 lobes desired beam pattern
p=sin(pi/180*-30);
ph=sin(pi/180*-20);
p2=sin(pi/180*20);
ph2=sin(pi/180*30);  
r=((1.*(t>p).*(t<ph)))+(1.*(t>p2).*(t<ph2));
figure, plot(tt,r)
%% 3 lobes desired beam pattern
% in this case put M1 = -50; M2 = 60;

% p=sin(pi/180*-50);
% ph=sin(pi/180*-40);
% p2=sin(pi/180*-15);
% ph2=sin(pi/180*15);  
% p3=sin(pi/180*40);
% ph3=sin(pi/180*60); 
% r=((1.*(t>p).*(t<ph)))+(1.*(t>p2).*(t<ph2))+(1.*(t>p3).*(t<ph3));
%% steering vector
A = zeros(L,N^2);
for l = 1:L
    aT = exp(  1j * pi*(0:N-1) * sin(theta(l)) ).'; % steering vector
    A(l,:) = kron(aT.',aT');
end
cvx_begin quiet
%%  Coefficient Assignment
 variable co(m)       
for i = 1:m
    su0(:,:) = pp0(i,:,:);  
    o = co(i) * su0 + o;   
end
tic
%% Optimization algorithm
minimize(norm(r - A * vec(o),2))
 subject to
    co>= 0;
cvx_end
toc;
%% Plot
P_sqp=zeros(1,L);
for l = 1:L
    et = exp(1j * pi * (0:N-1) * sin(theta(l))).';
    P_sqp(l) = et' * o * et;
end
plot(tt,10 * log10 (abs(P_sqp))), grid on, hold on
%% MSE
MSE = 1/L * ((r - P_sqp')' * (r - P_sqp'));
