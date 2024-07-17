% Stoica onProbing Min SLL

clc
clear all
% for counter = 1:100
delta = 1;
theta = (-90:delta:90)';
L = length(theta);
tic
%% Fig8a Stoica
% N=10;teta0=0;teta1=-10;teta2=10; %%SS(l)=A0*vec(R)-A(l,:)*vec(R)*abs(1e12*(omega(l)))^.1-t
% omega1=[-90:1:-20];omega2=[20:1:90];omega=[omega1,omega2];
%% Fig2c JianLi 
N = 16;
teta0 = 0;
teta1 = -6;
teta2 = 6;
omega1 = [-90:delta:-12];
omega2 = [12:delta:90];
omega = [omega1,omega2]; 
L1 = length(omega);

%% SQP Method
A = zeros(L1,N^2);
aT = exp(  1j * pi*(0:N-1)*sin((teta0*(pi/180)))).'; A0 = kron(aT.',aT');
aT = exp(  1j * pi*(0:N-1)*sin(teta1*(pi/180))).'; A1 = kron(aT.',aT');
aT = exp(  1j * pi*(0:N-1)*sin(teta2*(pi/180))).';A2 = kron(aT.',aT');
cvx_begin quiet
variable R(N,N) complex 
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

    %time (counter) = toc
    toc
    %clearvars -except time counter
% end
P_sqp = Beam_Pattern (N,theta*pi/180,R);
 %% 
plot(theta ,10*log10(abs(P_sqp)/max(P_sqp))),hold on, grid on 
