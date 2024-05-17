function [new_Volt]=inference_pred(Ind_Infer_xx, Ind_Infer_not, Bus_power_tk1, Bus_power_tk2, Volt_Bus_tk1, Part_Vi_Pj_xx, Part_Vi_Qj_xx)

%% Selecting part of the Z_matrix where row represents the the Infer_Buses and columns denotes the buses w.r.t. Inference is made....

Svp_xx = Part_Vi_Pj_xx(Ind_Infer_xx, Ind_Infer_not);
Svq_xx = Part_Vi_Qj_xx(Ind_Infer_xx, Ind_Infer_not);

%% Voltage Inference Multivariate Taylor series .....................

Jacob_delta = [ Svp_xx Svq_xx ];
   dp_xx = real(Bus_power_tk2) - real(Bus_power_tk1);
   dq_xx = imag(Bus_power_tk2) - imag(Bus_power_tk1);
    
delta_xx = [ dp_xx(Ind_Infer_not); dq_xx(Ind_Infer_not) ];

%% Updating voltages at the Infered buses ........................

new_Volt = Volt_Bus_tk1( Ind_Infer_xx ) + Jacob_delta*delta_xx;


end