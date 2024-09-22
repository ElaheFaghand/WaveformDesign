% The method of Basic Sets for symetric and asymetric
% without for cycles it works correctly 18/7/2024
% Revision has been strated in 18/9/2024
clc
clear 
tic
gamma = (-1:0.001:1);  % 2001 points are considered
theta = (asin(gamma)) ; 
% theta = (-90:.1:90) *pi/180;
K = length(theta);
M = 16;
w = ones(length(theta),1);
%% desired beampattern
start = -20 * pi/180; 
stop= 20 * pi/180;
pl1 = sin(start);
ph1 = sin( stop);
p_d = (1.*(sin(theta)>pl1).*(sin(theta)<ph1)); 
%% weight definition
% weight on sidelobes
for k = 1:K
    if theta(k) >= (stop + asin(1.5/M)) 
         w(k) = 1;
    end
    if theta(k) <= (start - asin(1.5/M))
        w(k) = 1;
    end
end
 %% 
A = zeros(K,M^2);
for l = 1:K
    a = exp( 1j * pi * (0:M-1) * gamma(l));
    A(l,:) = kron(a, conj(a));
end
%%
% J = zeros(M^2);
% for j = 1:M^2
%     for l = 1:M^2
%         for k = 1:K
%             J1 (k)  =  w(k) * A(k,j) * A(k,l)  ;
%         end
%         J (j,l) =  sum (J1);      
%     end
% end
% J  =  (A(:,1:M^2)).' * A(:,1:M^2) ;
Aw = A.*w;
J = Aw'*Aw;
%%
% for l = 1 : M^2
%     for k = 1:K
%         B1(k)  = w(k) * A(k,l)* p_d(k) ;
%     end
%     B(l) = sum(B1);
% end
b = (A'*(p_d'.*w));
 %%
%  r = inv(J) * B;
N = rank(J);
[V,D] = eig(J);
D = diag(D);
D = D(end:-1:length(D)-N+1); % for removing the zero ones
V = V(:,end:-1:length(V)-N+1); % for removing the zero ones
%%
r = V * ((V' * b)./D);
% r = zeros(M^2,1);
% for n = 1 : rank(J)
%     r = r + V(:,n) * (V(:,n)' * b)/(D(n,n));
% end

%%
RM = reshape(r, M, M);
% RM = RM+eye(M)*(max(max(abs(RM)))-RM(1));
P = Beam_Pattern (M, theta, RM);
% figure
% plot(theta*180/pi,10*log10(abs(P))), grid on, hold on
plot(theta*180/pi,(P)), grid on, hold on

eig(RM)
figure;
imagesc(abs(RM))
%%
% xline(start * 180/pi); xline(stop * 180/pi); yline(-3);
% MSE  = 1/K * ((p_d - P) * (p_d - P)');
toc