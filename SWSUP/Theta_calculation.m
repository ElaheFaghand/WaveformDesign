function  [eta, sum_Beta_i_P_i] = Theta_calculation (gamma, P_d , M, eta, beta)
    K = length(gamma);
    rc = .8;
%     w = 1; 
     %% 
%      for l = 1 : length(eta) 
         numerator = 1;  
             for w =1:100
%               while abs(numerator) > 1e-5

                  x = (gamma' -  eta);
                  P_shift = (sin (M*pi*(x)/2)./sin(pi*(x)/2)).^2/M^2;
                  P_shift(isnan(P_shift))= 1;

                  sum_Beta_i_P_i =  P_shift * beta;   
                   
                  rond1 = First_derivative (M, x);
                  rond2 = Second_derivative (M, x);

                  numerator =  (sum_Beta_i_P_i' - P_d) * rond1;
                  denumerator =  beta' .*( ones(1,K)* rond1.^2 ) +  ((sum_Beta_i_P_i' - P_d) * rond2);
                  eta = eta - rc * numerator./denumerator;
%                w = w+1;
              end
%      end
     %%
%     for j = 1 : length(eta) 
%         Rss(:,:,j) = shift_R_construction (eta(j), M, R_0);
%         P_R_ss(:,j) = Beam_Pattern (M, gamma,Rss(:,:,j));
%         plot(gamma , 10*log10(abs(P_R_ss(:,j)))), hold on, grid on, axis([-1, 1, -40,1]);
%         P_R_ss_weight(:,j) = beta(j) * P_R_ss(:,j); 
%     end
%     P_R_total = sum(P_R_ss_weight,2);
  sum_Beta_i_P_i =  P_shift * beta;   

%    plot(gamma , 10*log10(abs(sum_Beta_i_P_i))), hold on, grid on, axis([-1, 1, -40,1]);

end