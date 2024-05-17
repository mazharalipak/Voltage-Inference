%% Observability at the previous time step tk2 ...............

function [Volt_Bus_tk2, Bus_power_tk2] = obser_t2k(ii_indx_n, Ind_Infer_not, multi_att_P, multi_att_V, Ind_Infnot_att, tk2_a)

%          vecxx = ones(size(Ind_Infer_not,1),1);
%         vecxx([4,5,6],:) = multi_att;

%%%%%%%%%%%%%%%%%%%%% for creating time series plots %%%%%%%%%%%%%%%%%%%%%%

  Volt_Bus_tk2 = load('Volt_Bus_tk2');                                   % data of V(t_k+1) ..........................
  Volt_Bus_tk2 = Volt_Bus_tk2.Volt_Bus_tk2(:,tk2_a);                     % data of V(t_k+1) ..........................  
  Volt_Bus_tk2(Ind_Infnot_att)= multi_att_V.*Volt_Bus_tk2(Ind_Infnot_att);

        Power_tk2 = load('Power_tk2');                                   % data of S(t_k+1) ......................
        Power_tk2 = Power_tk2.Power_tk2(:,tk2_a);
    Bus_power_tk2 = Power_tk2(ii_indx_n);                                % S(t_k+1) corresponding to all the buses/present phases......
Bus_power_tk2(Ind_Infnot_att) = multi_att_P.*Bus_power_tk2(Ind_Infnot_att);
%%%%%%%%%%%%%%%%%%%%%% for creating surface %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Volt_Bus_tk2 = load('Bus_tk2_voltsss');
% Volt_Bus_tk2 = cell2mat( Volt_Bus_tk2.Bus_tk2_voltsss(index_Second) );
% Volt_Bus_tk2 = Volt_Bus_tk2(:,tk2_a);                                      % data of V(t_k+1) ..........................  
% Volt_Bus_tk2(Ind_Infer_not)= multi_att*Volt_Bus_tk2(Ind_Infer_not);
% 
% Power_tk2 = load('Power_tk2_voltsss');                                     % data of S(t_k+1) ......................
% Power_tk2 = cell2mat( Power_tk2.Power_tk2_voltsss(index_Second) );
% Power_tk2 = Power_tk2(:,tk2_a);
% Bus_power_tk2 = Power_tk2(ii_indx_n);                                      % S(t_k+1) corresponding to all the buses/present phases......

end