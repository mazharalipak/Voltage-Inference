%% Voltage Inference Framework ..................

clc;
clearvars;   
%% Optimization initialization .....................

alpha_par = 0.0005;                                                        % optimal step-size alpha ........
 max_Iter = 2000;                                                          % maximum iteration count .....
att_V_xx = ones(1000,1);                                                   % manipulating Vt_2 measurements at observable buses ...........
att_P_xx = ones(1000,1);

Observ_Buses = [ 4 10]; %[34 57];

%% 123 Bus  %[27 34 45 57]; %[3 14 68 72]; %[3 14 27 34 45 57 68 72]; %[10 20 60 80]; %[30 40 50 70]; 
 kappa_bus   = [2];

% attacked bus ................. 

 attack_bus = 5;                

%% choice of attack either voltages or nodal injections at the attack_bus .............

att_V_xx(100:100:1000) = ones(max(size(100:100:1000)),1); %+ (0.01:0.01:0.1)';
%0.25*(0.05:0.0055:0.1)';

att_P_xx(100:100:1000) = ones(max(size(100:100:1000)),1) + (0.01:0.01:0.1)';
% + 1*(0.05:0.0055:0.1)';

%% Initialization .................................
     mpc = load('IEEE_13.mat');                                            % data file IEEE test feeder..............
  Phases = phases_Ieee13;                                                  % phases/Bus ............except substation bus .......
           load('BIBC.mat')
           load('BCBV.mat')
           load('Part_Vi_Pj_xx');                                          % Senstivities of V_i w.r.t. P_j .........
           load('Part_Vi_Qj_xx');                                          % Senstivities of V_i w.r.t. Q_j .........

  n_branch = size(struct2table(mpc.Branch),1);
[n_branch, Bus_no_xx, ii_indx_n, zz_mat, Vs, Bus_phases,...                % Initalization data file ..............
    phases_xx] = inference_data(n_branch, Phases);

%% Z = BCBV(BIBC) matrix .......................

%  [C_mat, D_mat, E_mat] = svd(zz_mat.zz_mat);
%              yy_mat_xx = E_mat*inv(D_mat)*C_mat';

               yy_mat_xx = inv(zz_mat.zz_mat);          
%% Buses where Inference is made ...... locations of Infer_Buses with their corresponding phases .........................................

% buses where voltages are infered ................  

 Infer_Buses = 2:14;
Infer_Buses(Observ_Buses-1) = []; 

   Infer_indx = ismember( Bus_phases, Infer_Buses );                       % locations where Bus_phases equivalent to Infer_Buses ........
 Ind_Infer_xx = find(Infer_indx);                                          % locations of buses/phases with unobservability and requires Inference........
Ind_Infer_not = find(~Infer_indx);                                         % locations of buses/phases with observability........

Ind_Infnot_att = find(ismember( Bus_phases, attack_bus ));                 % Finding location of the attacked bus ..........                        

% Location of the overlapping bus/phases ...........................

   find_kappa = find(ismember(Bus_phases, kappa_bus));
 Ind_Infer_xx = find_kappa;                                                % inferring only at the kappa locations ....
 
%% Observability at the previous time step tk1 for all buses ...............

[Volt_Bus_tk1, Bus_power_tk1] = obser_t1k(ii_indx_n);                      % for all buses with corresponding phases except slack ........

  
 for i_inner = 1:1000                                                       % Iteration counter for each cell from the V(t_k) and S(t_k) from the obser_tk2 script ........
 
tk2_a = i_inner;                                                          % Internal Iteration counter  for each cell ...........

 %%

 multi_att_V = att_V_xx(i_inner)*ones( (max(size(Ind_Infnot_att))),1 );
 multi_att_P = att_P_xx(i_inner)*ones( (max(size(Ind_Infnot_att))),1 );

%% Observability at the next time step tk2 .............................

[Volt_Bus_tk2, Bus_power_tk2] = obser_t2k(ii_indx_n, Ind_Infer_not,...
    multi_att_P, multi_att_V, Ind_Infnot_att, tk2_a);                                       % for all buses with corresponding phases except slack ........

%Bus_power_tk2( [20, 21, 22] ) = Bus_power_tk2( [20, 21, 22] );            % altering measurements ...................

%% Prediction step ................

[new_Volt] = inference_pred(Ind_Infer_xx, Ind_Infer_not, Bus_power_tk1,... % multivariate Taylor series prediction step ....... 
    Bus_power_tk2, Volt_Bus_tk1, Part_Vi_Pj_xx, Part_Vi_Qj_xx); 

[error_pred1, error_pred2] = pred_plots(Ind_Infer_xx, Volt_Bus_tk2,...     % prediction plots & errors .........
    new_Volt);

%% Corrector step estimating complex voltages ...................

               Vs_vecx = cell(n_branch,1);
Vs_vecx (1:n_branch,1) = {Vs};                                                                               
                Vs_vec = cell2mat(Vs_vecx);

%% Optimization ..........................

            Updated_Volt = Volt_Bus_tk2;

%% Fixed Point Iterations code .....................................

Tol   = 1;  
Iter  = 1;
count = 0;

while (Tol > 1e-12)

new_Volt1 = new_Volt;

Updated_Volt(Ind_Infer_xx) = new_Volt;

                  Delta_Vx = Vs_vec(ii_indx_n) -  Updated_Volt;
                   matxx_V = diag( conj(Updated_Volt) ) ;
                estimate_P =  ( matxx_V*yy_mat_xx )*Delta_Vx ;
                     P_xxc =  estimate_P(Ind_Infer_not)...
                           - conj( Bus_power_tk2(Ind_Infer_not) );

% Objective function 1/2|| S_b*( Va(t_k+1) ) - S_b*(t_k+1) || ..........

error_P = 0.5*([ real(P_xxc); imag(P_xxc) ]'*[ real(P_xxc); imag(P_xxc) ]);

%Grad_x_new = -( matxx_V(Ind_Infer_not, Ind_Infer_not)*yy_mat_xx(Ind_Infer_not, Ind_Infer_xx) )'*P_xxc;
   
%% Calculating gradient ...........................

test_xx1 = diag(yy_mat_xx*Delta_Vx) - ( conj(Updated_Volt).*yy_mat_xx );               % partial Si / partial V^r_j
test_xx2 = -1i*( diag(yy_mat_xx*Delta_Vx) + ( conj(Updated_Volt).*yy_mat_xx ) );       % partial Si / partial V^m_j

% Jacob = [ (P_i/partial V^r_j) (P_i/partial V^m_j); (Q_i/partial V^r_j) (Q_i/partial V^m_j) ] ...............

Jacob_deltaxx = [ real(test_xx1(Ind_Infer_not, Ind_Infer_xx))...
    real(test_xx2(Ind_Infer_not, Ind_Infer_xx));...
    imag(test_xx1(Ind_Infer_not, Ind_Infer_xx))...
    imag(test_xx2(Ind_Infer_not, Ind_Infer_xx)) ];

  Grad_x_new  = Jacob_deltaxx'*[ real(P_xxc); imag(P_xxc) ];
%% Optimal Step-Size Armijo line search ...............

[alphas] = alpha_grad(error_P, alpha_par, ii_indx_n, Ind_Infer_xx,...
    Ind_Infer_not, Volt_Bus_tk2, new_Volt, yy_mat_xx, Grad_x_new, Vs_vec, Bus_power_tk2);

%% Updaing the vector 

delta_xv =  alpha_par*Grad_x_new;
new_Volt = new_Volt - ( delta_xv( 1: length(delta_xv)/2 ) +...
    1i*( delta_xv( (length(delta_xv)/2)+1: end ) ) );

%%  Gradient Decent stoping criteria ...................

error_volt = [ ( abs(real(new_Volt)) - abs(real(new_Volt1)) );...
    ( abs(imag(new_Volt)) - abs(imag(new_Volt1)) )];

%% Storing different errors during gradient descent iterations ..............

      loss_f(Iter) = error_P;                           
 grad_error (Iter) = max( abs( [real(Grad_x_new); imag(Grad_x_new)] ) );
Volt_errors (Iter) = max( abs( error_volt ) );
 Iter_alpha (Iter) = alphas;

%% Algorithm closing structure ............................

Iter = Iter + 1;
Tol =  max( abs( error_volt ) );                    
count = count + 1;
    if count == max_Iter
        break;
    end
end


%% Using Inference recursively for the inner iterations .... loop t_{k-1} used for t_k .........

               Volt_Bus_tk1 = Volt_Bus_tk2;
 Volt_Bus_tk1(Ind_Infer_xx) = new_Volt;

             Bus_power_tk1 = Bus_power_tk2;
Bus_power_tk1(Ind_Infer_xx)= conj( estimate_P(Ind_Infer_xx) );

% Power...........

  %            Bus_voltxx_tk1 = zeros(3*n_branch, 1);
  % Bus_voltxx_tk1(ii_indx_n) = Volt_Bus_tk1;
  % 
  % Bus_powerith_tk1 = zeros(3*n_branch, 1);
  % Bus_powerith_tk1(ii_indx_n) = Bus_power_tk1;
  % 
  % [ J_full, J11, J12 ] = tnr_Jacob( ii_indx_n, BCBV, BIBC, Bus_powerith_tk1, Bus_voltxx_tk1 );
  % 
  % [Part_Vi_Pj_xx, Part_Vi_Qj_xx, magV_Pj_xx, magV_Qj_xx] = senstivity_Dist(zz_mat.zz_mat,...
  %    ii_indx_n, Bus_voltxx_tk1, J11, J12 );

%% Plots for correction step  ...............

[error_corr1, error_corr2, Table_Inf]=corr_plots(Ind_Infer_xx, Bus_phases,...
    phases_xx, Volt_Bus_tk2, new_Volt);

%% Calculating errors at each inner iterations t_k for P_b and Q_b ..............

          PP_infer_xx = conj( estimate_P );
Pinf_xx = max(abs( abs( real( PP_infer_xx(Ind_Infer_xx) ) ) -...
    abs( real( Bus_power_tk2(Ind_Infer_xx) ) ) ) );

Qinf_xx = max(abs( abs( imag( PP_infer_xx(Ind_Infer_xx) ) ) -...
    abs( imag( Bus_power_tk2(Ind_Infer_xx) ) ) ) );

%% Storing data at each inner iteration .....................

error_corr1_tk2(i_inner) = error_corr1;                                    % Storing absolute difference of ( max| |\hat{V}_b (t_k)|    - |V_b(t_k)| | ) 
error_corr2_tk2(i_inner) = error_corr2;                                    % Storing absolute difference of ( max| |\hat{\delta}_b (t_k)|    - |\delta_b(t_k)| | ) 

   error_P1_tk2(i_inner) = Pinf_xx;                                        % Storing absolute difference of ( max| |\hat{P}_b (t_k)|    - |P_b(t_k)| | ) 
   error_Q1_tk2(i_inner) = Qinf_xx;                                        % Storing absolute difference of ( max| |\hat{Q}_b (t_k)|    - |Q_b(t_k)| | ) 

P_Infer_act(:,i_inner) = abs( real(Bus_power_tk2(Ind_Infer_xx)) ); 
P_Infer_est(:,i_inner) = abs( real(PP_infer_xx(Ind_Infer_xx)) ); 

Iterr_inference(i_inner) = Iter;
  Tol_inference(i_inner) = Tol;

V_comp_kappa(:,i_inner) = new_Volt;                                        % Complex voltage phasor at unobservable location \hat{V}(t_k).........

V_kappa_obser(:,i_inner) = Volt_Bus_tk2(Ind_Infnot_att);                     % Complex voltage phasor at observable location V(t_k).........

%(find_kappa);        

%max(error_corr1_tk2)
max(error_corr2_tk2)

max(Iterr_inference)

i_inner

 end

