%% Stoica 2007 SQP Method for solving SLL design problem in section III-D for wide band beampatterns in article: 
%% On Probing Signal Design For MIMO Radar
%% It is written by Elahe Faghand
clc
clear all
tic
theta = (-90:1:90)'*(pi/180);
L = length(theta);
N = 16; % Number of antennas 
%% desired Beam pattern angles
teta0 = 0; 
teta1 = -20 - 3;
teta2 = 20 + 3;
teta = [teta1+3:1:teta2-3];
omega1 = [-90:1:-28];
omega2 = [28:1:90];
omega = [omega1,omega2];
L1 = length(omega);

%% SQP Method
% steering vector definitons
A = zeros(L1,N^2);
aT = exp( 1j * pi*(0:N-1)*sin(teta0*(pi/180))).';
A00 = kron(aT.',aT');
aT = exp( 1j * pi*(0:N-1)*sin(teta1*(pi/180))).';
A1 = kron(aT.',aT');
aT = exp( 1j * pi*(0:N-1)*sin(teta2*(pi/180))).';
A2 = kron(aT.',aT');

cvx_begin quiet
variable R(N,N) complex  % covariance matrix
variable t 
minimize(-t)
subject to
  R == hermitian_semidefinite(N);
  R(logical(eye(N))) == 1/N;
for k=1:length (teta)
    aT = exp(  1j * pi*(0:N-1)*sin(teta(k)*(pi/180))).';
    A0(k,:) = kron(aT.',aT');
    for l = 1:L1
        aT = exp(  1j * pi*(0:N-1)*sin(omega(l)*pi/180) ).';
        A(l,:) = kron(aT.',aT');
        SS(l) = A0(k,:) * vec(R) - A(l,:)*vec(R) - t;
    end
    real(SS) >= 0;
    k;
end
A1*vec(R)-0.5*A00*vec(R) == 0;  % for multi lobes i removed this line
A2*vec(R)-0.5*A00*vec(R) == 0;  
cvx_end
toc
R =  R/sum(sum(R));
P_sqp = Beam_Pattern(N, theta, R); % achieved beampattern
plot(theta*180/pi ,10*log10(abs(P_sqp))),grid on,hold on
%%
function Beam_Pattern = Beam_Pattern (N,theta,R)
Beam_Pattern = zeros(1, length(theta));
    for l = 1:length(theta)
        et = exp(1j*pi*(0:N-1)*sin(theta(l))).';
        Beam_Pattern(l) = et'*R*et;
    end
end
