% Stoica onProbing Min SLL
clc
clear all
tic
theta = (-90:1:90)'*(pi/180);
L = length(theta);

%% Fig 3a JianLi
N = 16;
teta0 = 0; 
teta1 = -20 - 3;
teta2 = 20 + 3;
teta = [teta1+3:1:teta2-3];
% multi lobe
% teta1 = 35 - 3;
% teta2 = 55 + 3;
% teta_2 = [teta1+3:1:teta2-3];
% 
% teta1 = -55 - 3;
% teta2 = -35 + 3;
% teta_3 = [teta1+3:1:teta2-3];
%  
% teta = [teta teta_2 teta_3];
%%
omega1 = [-90:1:-28];
omega2 = [28:1:90];
omega = [omega1,omega2];

% multi lobe

% omega1 = [-90:1:-61];
% omega2 = [-29:1:-16];
% omega3 = [16:1:29];
% omega4 = [61:1:90];
% % 
% omega = [omega1 omega2 omega3 omega4];
L1 = length(omega);

%% SQP Method
A = zeros(L1,N^2);
aT = exp( 1j * pi*(0:N-1)*sin(teta0*(pi/180))).';
A00 = kron(aT.',aT');
aT = exp( 1j * pi*(0:N-1)*sin(teta1*(pi/180))).';
A1 = kron(aT.',aT');
aT = exp( 1j * pi*(0:N-1)*sin(teta2*(pi/180))).';
A2 = kron(aT.',aT');

cvx_begin quiet
variable R(N,N) complex 
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
        % SS(l)=A0(k,:) *vec(R)-A(l,:)*vec(R)*abs(2e2*(omega(l)))^0.005-t;
    end
    real(SS) >= 0;
    k;
end
A1*vec(R)-0.5*A00*vec(R) == 0;  % for multi lobes i removed this line
A2*vec(R)-0.5*A00*vec(R) == 0;  
cvx_end
toc
R=  R/sum(sum(R));
P_sqp = Beam_Pattern(N, theta, R);
plot(theta*180/pi ,10*log10(abs(P_sqp))),grid on,hold on
