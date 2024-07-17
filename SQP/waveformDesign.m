%% wave form
clc
clear
theta = (-90:0.1:90)' *pi/180;
tt = (-90:0.1:90);
N = 10;
L = 1024;
% Rg = ones(N);
Rg = importdata("R.mat");
P = Beam_Pattern (N,theta,Rg);
plot(tt,10*log10(abs(P)/max(abs(P)))); grid on, hold on,
%%
[u,v] = eig(Rg);
si = randn(L,N);
X = (si* v^0.5 *u');
% R_BPSK=(2/pi)*asin(Rg);
% R = sin((pi/2)*Rg);
% eig(R);
% ZBPSK=sign(X)/sqrt(2);
Z_BPSK = sign(X)/sqrt(2);
X_BPSK = (Z_BPSK'*Z_BPSK);
R_BPSK = X_BPSK/L;
P_BPSK = Beam_Pattern (N,theta,R_BPSK);
plot(tt,10*log10(abs(P_BPSK)/max(abs(P_BPSK)))); grid on, hold on,
%% QPSK
X = reshape(X,1,N*L);
for i = 1:N*L/2
    ZQPSK(i) = (sign(X(2*i-1)) + 1j*sign(X(2*i)))/sqrt(2);
end
Z_QPSK = reshape(ZQPSK,L/2,N);
%  stem(ZQPSK,'o');
%  stem(ZBPSK,'*')
X_QPSK = (Z_QPSK'*Z_QPSK);
R_QPSK = X_QPSK/L;
P_QPSK = Beam_Pattern (N,theta,R_QPSK);
plot(tt,10*log10(abs(P_QPSK)/max(abs(P_QPSK)))); grid on, hold on,
%% satrha:symbol    sotton:anten
% autocorr(ZBPSK(:,1));  % same below, positive side
crosscorr(Z_BPSK(:,1),Z_BPSK(:,1));   % anten 1 
crosscorr(Z_BPSK(:,1),Z_BPSK(:,3));   % anten 1  & anten2
crosscorr( real(Z_QPSK(:,1)), real(Z_QPSK(:,1))); % anten 1 
crosscorr( real(Z_QPSK(:,1)), real(Z_QPSK(:,2))); % anten 1  & anten2
%%

%%
function Beam_Pattern = Beam_Pattern (N,theta,R)
Beam_Pattern = zeros(1, length(theta));
    for l = 1:length(theta)
        et = exp(1j*pi*(0:N-1)*sin(theta(l))).';
        Beam_Pattern(l) = et'*R*et;
    end
end
