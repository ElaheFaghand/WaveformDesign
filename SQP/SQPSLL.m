%% Stoica 2007 SQP Method for solving SLL design problem in section III-D in article: 
%% On Probing Signal Design For MIMO Radar
%% It is written by Elahe Faghand

clc
clear all
delta = 1;
theta = (-90:delta:90)';
tic
%% Fig8a Stoica properties
N=10; % number of antennas
teta0=0; % centeral angle
teta1=-10;
teta2=10; 
omega1=[-90:1:-20];
omega2=[20:1:90];
omega=[omega1,omega2];
L1 = length(omega);
%% SQP Method
A = zeros(L1,N^2);
aT = exp(  1j * pi*(0:N-1)*sin((teta0*(pi/180)))).'; A0 = kron(aT.',aT'); % steering vector
aT = exp(  1j * pi*(0:N-1)*sin(teta1*(pi/180))).'; A1 = kron(aT.',aT');
aT = exp(  1j * pi*(0:N-1)*sin(teta2*(pi/180))).';A2 = kron(aT.',aT');
cvx_begin quiet
variable R(N,N) complex  % covariance matrix
variable t 
minimize(-t)
subject to
  R == hermitian_semidefinite(N);
  R(logical(eye(N))) == 1/N;
for l = 1:L1
    aT = exp(  1j * pi*(0:N-1)*sin(omega(l)*pi/180) ).';
    A(l,:) = kron(aT.',aT');
    SS(l) = A0*vec(R) - A(l,:)*vec(R) - t;
   % SS(l)=A0*vec(R)-A(l,:)*vec(R)*(5000*weight(l))-t;
end
real(SS)>=0;
A1*vec(R) - 0.5*A0*vec(R) == 0;
A2*vec(R) - 0.5*A0*vec(R) == 0;  
cvx_end
toc
P_sqp = Beam_Pattern (N,theta*pi/180,R);
plot(theta ,10*log10(abs(P_sqp)/max(P_sqp))),hold on, grid on 
%%
function Beam_Pattern = Beam_Pattern (N,theta,R)
Beam_Pattern = zeros(1, length(theta));
    for l = 1:length(theta)
        et = exp(1j*pi*(0:N-1)*sin(theta(l))).';
        Beam_Pattern(l) = et'*R*et;
    end
end
