% Shift Evaluation
function  [eta, sum_Beta_i_P_i] = Theta_calculation (gamma, P_d , M, eta, beta)
    K = length(gamma);
    rc = .8;
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
  sum_Beta_i_P_i =  P_shift * beta;   
end
