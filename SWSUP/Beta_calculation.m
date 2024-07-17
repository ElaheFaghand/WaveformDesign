function  [p_L_Final, beta, P_shift]  = Beta_calculation (gamma, P_d, M, eta, beta_initial)
    K = length(gamma);
    I =  length(eta);
    P_shift = zeros(K,I);

    %% formul
    x = (gamma' - eta);
    P_shift = (sin (M*pi*x/2)./sin(pi*x/2)).^2/M^2;
%     for k=1:K
%         for l =1:I
%             P_shift (k,l) = (sin (M*pi*(gamma(k) - eta(l))/2)/sin(pi*(gamma(k) - eta(l))/2))^2/M^2;
%         end
%     end
    P_shift(isnan(P_shift))= 1;
    p_L =  P_shift * beta_initial';   
            %%
%     for i = 1:I
%           R_shift(:,:,i) = shift_R_construction (eta(i) , M , R_0);
%           P_shift(:,i) = Beam_Pattern (M , gamma , R_shift(:,:,i)); 
%           plot(gamma , 10*log10(abs(P_shift(:,i)))), hold on, grid on, axis([-1, 1, -40,1]);
%           % plot(theta , 10*log10(abs(P_shift(:,i)))), hold on, grid on, axis([-90, 90, -40,1]);
%           P_shift_weight(:,i) = beta_initial(i) * P_shift(:,i); 
%     end
% p_L = sum(P_shift_weight,2);
    %% Gradian
%     for l = 1:I
%             for k = 1 : K
%                 G1(k) = -2 * (P_d(k) - sum_P_shift_weight(k)) * P_shift(k,l);            
%             end
%         G(l) = sum (G1);
%     end

 G = -2 * (P_d' - p_L)' * P_shift;            
 
    %% Hessian Matrix calculation
%     H = zeros(I);
%      for l = 1:I
%          for j = 1:I
%              for k = 1:K
%                  H1 (l,j,k) = 2 * P_shift(k,l) * P_shift(k,j); 
%              end
%          end
%      end
%      H = sum (H1 , 3);
    H = 2 * P_shift' * P_shift;
    beta = beta_initial' - inv((H)) * G' ;
    %% plot
    p_L_Final =  P_shift * beta;   
%   plot(theta, 10*log10(abs(p_L_Final))), grid on, hold on, axis([-80 80, -40 ,1]); 
end