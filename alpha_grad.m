%% Armijo inexact line search
function [alphas] = alpha_grad(pfxx_vec, alpha_par, ii_indx_n, Ind_Infer_xx, Ind_Infer_not, Volt_Bus_tk2, new_Volt, yy_mat_xx, Grad_x_new, Vs_vec, Bus_power_tk2)

% Descent direction - gradient fun ...............

grad_fun = - ( Grad_x_new( 1: (length(Grad_x_new)/2) ) +...
     1i* Grad_x_new( (length(Grad_x_new)/2) + 1: end ) );

% Finding f(x + alpha*grad_fun) ...............
               Updated_Volt = Volt_Bus_tk2;
                 new_vector = new_Volt + alpha_par*grad_fun;
Updated_Volt (Ind_Infer_xx) = new_vector;
                    Delta_Vx = Vs_vec(ii_indx_n) -  Updated_Volt;
                     matxx_V = diag( conj(Updated_Volt) ) ;
                  estimate_P =  ( matxx_V*yy_mat_xx )*Delta_Vx ;
                       P_xxc =  estimate_P(Ind_Infer_not)...
                           - conj( Bus_power_tk2(Ind_Infer_not) );
                    P_xxc_xx =  0.5*( [real(P_xxc);imag(P_xxc)]'*[real(P_xxc);imag(P_xxc)]);

if  P_xxc_xx < ( pfxx_vec + (10^-4)*alpha_par*(grad_fun'*grad_fun) ) 
    alphas = alpha_par;
else
    alphas = alpha_par/2;
end

end
