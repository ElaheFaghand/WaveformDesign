% The method of Basic Sets for symetric and asymetric
% without for cycles it works correctly 18/7/2024
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
start = -10 * pi/180;
stop= 30 * pi/180;
pl1 = sin(start);
ph1 = sin( stop);
p_d = (1.*(sin(theta)>pl1).*(sin(theta)<ph1)); 
%
% start_1 =  -58*pi/180;
% stop_1 = -32*pi/180;
% 
% start_2 = -13*pi/180;
% stop_2 = 13*pi/180;
% 
% start_3 = 32*pi/180;
% stop_3 = 58*pi/180;
% 
% p = sin(start_1);
% ph = sin(stop_1);
% 
% p2 = sin(start_2);
% ph2 = sin(stop_2);  
% 
% p3 = sin(start_3);
% ph3 = sin(stop_3);  
% t = sin(theta);
% p_d = ((1.*(t>p).*(t<ph)))+(1.*(t>p2).*(t<ph2))+(1.*(t>p3).*(t<ph3));
%% weight definition
% weight on sidelobes
for k = 1:K
    if theta(k) >= (stop + asin(1.5/M)) 
         w(k) = 1000;
    end
    if theta(k) <= (start - asin(1.5/M))
        w(k) = 1000;
    end
end
% 
% figure
% plot(theta*180/pi, w)

% weight on mainlobe
% for k = 1:K
%     if (theta(k)<=stop + 2/M) && (start - 2/M<= theta(k)) 
%           w(k) = 10000;
%     end
% end

% figure
% plot(theta*180/pi,w)
% 2 mainlobes
% for k = 1:K
%     if theta(k) >= (stop_1 + 1.5/M) && theta(k) <= (start_2 - 1.5/M) || theta(k) >= (stop_2 + 1.5/M) && theta(k) <= (start_3 - 1.5/M)
%          w(k) = 10;
%     end
% 
%     if theta(k) <= (start_1 - 1.5/M) %&& theta(k)>=(start - 1)
%         w(k) = 10;
%     end
% 
%     if theta(k) >= (stop_3 + 1.5/M) %&& theta(k)<=(stop + 5)
%          w(k) = 10;
%     end
% 
% end
 %% 
A = zeros(K,M^2);
for l = 1:K
    a = exp( 1j * pi * (0:M-1) * gamma(l));
    A(l,:) = kron(a, conj(a));
end
% A = importdata("A.mat"); 

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
%%

% J = zeros(M^2);
% for k = 1:K
%     J  =  J + (w(k) * A(k,:))' * w(k)*A(k,:)  ;
% end
% J1  =  ((conj(Aw)).')*(Aw);

Aw = A.*w;
Awh=(Aw)';
J  =  Awh * Aw;

%%
% for l = 1 : M^2
%     for k = 1:K
%         B1(k)  = w(k) * A(k,l)* p_d(k) ;
%     end
%     B(l) = sum(B1);
% end

% B =  (p_d * A(:,(1:M^2))); 

% b = zeros(M^2 , 1);
% for k = 1:K
%     % B =  B + w(k) * A(k,(1:M^2)).'* p_d(k) ;
%     b =  b + w(k) .* A(k,:)'.* p_d(1,k) ;
% end
b = (A'*(p_d'.*w));

 %%
 % R = inv(J) * B;
 % u = importdata("u.mat");
 % v = importdata("v.mat");
[V, D] = eig(J);
D = diag(D);
N = rank(J);
D = D(1:N); 
V = V(:, 1:N);
r = V * ((V' * b)./D);

% r = zeros(M^2,1);
% for n = 1 : rank(J)
%     r = r + V(:,n) * (V(:,n)' * B)/(D(n,n));
% end

% for k=1:2*M-1
%     r = r+ u(:,k)'* B/v(k,k)*u(:,k);
% end

%% it is wrong and more time consuming
% v = diag(v);
% alpha = (u(:,1 : 2*M)' * B')./v(1 : 2*M);
% r =  u(:,1 : 2*M) * alpha(1 : 2*M);
%%
RM = reshape(r, M, M);
RM = RM+eye(M)*(max(max(abs(RM)))-RM(1));
P = Beam_Pattern (M, theta, RM);
% figure
% plot(theta *180/pi,10*log(abs(P))), grid on, hold on
plot(theta*180/pi,10*log10(abs(P))), grid on, hold on

% figure;
% imagesc(abs(RM))

%% making a shift
% RM_shift = shift_R_construction(-10 , M, RM);
% P_shift = Beam_Pattern (M, theta, RM_shift);
% % figure
% plot(theta *180/pi, 10*log10(abs(P_shift))), grid on, hold on

%%
% xline(start * 180/pi); xline(stop * 180/pi); yline(-3);
% MSE  = 1/K * ((p_d - P) * (p_d - P)');
toc