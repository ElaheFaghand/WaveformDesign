%% UPA-LIB PAPER1
clc
clear
for ii = 1:100

    format long
    N = 30; % number of antennas
    of = 1; 
    i = 1;
%     M1 = -50; % strart point
%     M2 = 60;  % stop point
    theta_start1 = randi([0,40]);
    theta_stop1 = randi ([theta_start1+5, 80]);

    theta_stop = randi([-40,0]);
    theta_start = randi ([-80,theta_stop-5]);

    M1 = theta_start;
    M2 = theta_stop1;
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
       
    % for k = 0 : N-1
    %     if abs(p)>=0 
    %     o(k+1) = exp(1i * pi * k * tan(p * pi/180));
    %     else
    %     o(k+1)=0;
    % end
    % end
    
    %% plot members of UPA library
    % for l = 1:length (theta)
    %     et = exp(1j * pi * (0:N-1) * sin(theta(l))).';
    %     P_sqp(l) = et'* (o'*o) * et;
    % end
    pp0(i,:,:) = o'*o;
    i = i + 1;
    end
    
    m = m2 - m1 + 1;
    t = sin(theta);
    o = zeros(N);
    
%% 1- desired beam pattern
%     pl1 = sin(pi/180 * M1);
%     ph1 = sin(pi/180 * M2);
%     r = (1.*(t>pl1).*(t<ph1));
%% 2- desired beam pattern
p=sin(pi/180*theta_start1);
ph=sin(pi/180*theta_stop1);
p2=sin(pi/180*theta_start);
ph2=sin(pi/180*theta_stop);  
r=((1.*(t>p).*(t<ph)))+(1.*(t>p2).*(t<ph2));
% figure, plot(tt,r)
%% 3- desired beam pattern
%     p=sin(pi/180*-50);
%     ph=sin(pi/180*-40);
%     p2=sin(pi/180*-15);
%     ph2=sin(pi/180*15);  
%     p3=sin(pi/180*40);
%     ph3=sin(pi/180*60); 
%     r=((1.*(t>p).*(t<ph)))+(1.*(t>p2).*(t<ph2))+(1.*(t>p3).*(t<ph3));
    %% steering vector
    A = zeros(L,N^2);
    for l = 1:L
        aT = exp(  1j * pi*(0:N-1) * sin(theta(l)) ).';
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
    time(ii) = toc;

    %% Plot
    P_sqp=zeros(1,L);
    for l = 1:L
        et = exp(1j * pi * (0:N-1) * sin(theta(l))).';
        P_sqp(l) = et' * o * et;
    end
    % FigH = figure(1);
    % set(FigH,'visible',"off")
    % set(gcf,"Position",[1 1 2 1]);
    % 
    % set(FigH,'Units','Inches');
    % pos = get(FigH,'Position');
    % set(FigH,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%     figure
%     plot(tt,10 * log10 (abs(P_sqp))), grid on, hold on
    % saveas (FigH,'D:\mimo_papers_99_09_04\mimo_papers_99_09_04\matlab code\version14\UQP\UQP',"pdf")
    
    %% MSE
    MSE(ii) = 1/L * ((r - P_sqp')' * (r - P_sqp'));
%     ii= ii+1;
    clearvars -except ii MSE time
end
mean(MSE)
var (MSE)
time = time/1.17;
mean(time)
var (time)
% mean(time)