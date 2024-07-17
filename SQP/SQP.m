%% Stoica 2007 SQP Method

clc
clear all
% for i = 1 : 100
    format long
    
    % gamma = (-1:0.001:1);  % 2001 points are considered
    theta = (-90:0.1:90)' *pi/180;
    % theta = (asin(gamma))' ;
    tt = theta * 180/pi;
    L = length(theta);
    t=sin(theta);
    N = 10;
    %% 
    % p=sin(pi/180*-50);
    % ph=sin(pi/180*-40);
    % p2=sin(pi/180*-15);
    % ph2=sin(pi/180*15);  
    % p3=sin(pi/180*40);
    % ph3=sin(pi/180*60);
    % r=((1.*(t>p).*(t<ph)))+(1.*(t>p2).*(t<ph2))+(1.*(t>p3).*(t<ph3));
    %%
    pl1=sin(pi/180*-30);
    ph1=sin(pi/180*30);
    r=(1.*(t>pl1).*(t<ph1));
    %%
    % weight = ones(L,1);
    % weight (find(tt==-40):find(tt==-24)) =10;
    % weight (find(tt==24):find(tt==40)) =10;
    %%
    % p=sin(pi/180*-30);
    % ph=sin(pi/180*-10);
    % p2=sin(pi/180*10);
    % ph2=sin(pi/180*30);  
    % r=((1.*(t>p).*(t<ph)))+(1.*(t>p2).*(t<ph2));
    %% 
    % p=sin(pi/180*-70);
    % ph=sin(pi/180*-60);
    % p2=sin(pi/180*-40);
    % ph2=sin(pi/180*-30);  
    % p3=sin(pi/180*-10);
    % ph3=sin(pi/180*5); 
    % p4=sin(pi/180*30);
    % ph4=sin(pi/180*55); 
    % r=((1.*(t>p).*(t<ph)))+(1.*(t>p2).*(t<ph2))+(1.*(t>p3).*(t<ph3))+(1.*(t>p4).*(t<ph4));
    %%
    % p=sin(pi/180*-40);
    % ph=sin(pi/180*-35);
    % p2=sin(pi/180*-25);
    % ph2=sin(pi/180*-15);  
    % p3=sin(pi/180*-5);
    % ph3=sin(pi/180*5); 
    % p4=sin(pi/180*15);
    % ph4=sin(pi/180*20); 
    % p5=sin(pi/180*40);
    % ph5=sin(pi/180*50); 
    % r=((1.*(t>p).*(t<ph)))+(1.*(t>p2).*(t<ph2))+(1.*(t>p3).*(t<ph3))+(1.*(t>p4).*(t<ph4))+(1.*(t>p5).*(t<ph5));
    %% SQP Method
    A = zeros(L,N^2);
    for l = 1:L
        aT = exp(  1j * pi*(0:N-1)*sin(theta(l)) ).';
        A(l,:) = kron(aT.',aT');
    end
    
    cvx_begin quiet
    variable R(N,N) complex
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
     plot(tt,10*log10(abs(P_sqp)/max(abs(P_sqp))))
    %% desired plot
    % for i=1:length(r)
    %     if r(i)==0;
    %         r(i)=1e-2;
    %     end
    % end
    % plot(tt,10*log10(r))
    % time (i) = toc;
    % clearvars -except time i
% end
 % MSE  = 1/L * ((r' - P_sqp) * (r' - P_sqp)')
 
% mean(time)