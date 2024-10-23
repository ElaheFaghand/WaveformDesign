% Weight Evaluation
function  [p_L_Final, beta, P_shift]  = Beta_calculation (gamma, P_d, M, eta, beta_initial)
    K = length(gamma);
    I =  length(eta);
    P_shift = zeros(K,I);

    %% formul
    x = (gamma' - eta);
    P_shift = (sin (M*pi*x/2)./sin(pi*x/2)).^2/M^2;

    P_shift(isnan(P_shift))= 1;
    p_L =  P_shift * beta_initial';   
    %% Gradian
    G = -2 * (P_d' - p_L)' * P_shift;            
    %% Hessian Matrix calculation
    H = 2 * P_shift' * P_shift;
    beta = beta_initial' - inv((H)) * G' ;
    %% plot
    p_L_Final =  P_shift * beta;   
end
