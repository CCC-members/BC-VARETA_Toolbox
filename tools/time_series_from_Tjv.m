function [J_time] = time_series_from_Tjv(Tjv,V_time,IsCurv,IsField,scaleKe,GridOrient,GridAtlas)
%% obtain source time series (J_time) from data (V_time) using transfer operator Tjv
%% V_time was obtained from the hilbert envelope of the filtered MEG/EEG signal used to compute Tjv
Ng     = size(Tjv,1);
Nt     = size(V_time,2);
if IsCurv == 0
    J_time               = Tjv*V_time;
    if IsField == 2 || IsField == 3
        J_time           = reshape(J_time,3,Ng/3,Nt);
        J_time           = sqrt(squeeze(sum(abs(J_time).^2,1)));
        Tjv              = bst_gain_orient(Tjv',GridOrient,GridAtlas);
        J_time_angle     = exp(i*angle(Tjv*V_time));
        J_time           = J_time.*J_time_angle;
    end
elseif IsCurv == 1
    J_time               = squeeze(Tjv(:,:,1))*V_time + squeeze(Tjv(:,:,2))*V_time;
    if IsField == 2 || IsField == 3
        J_time           = reshape(J_time,3,Ng/3,Nt);
        J_time           = sqrt(squeeze(sum(abs(J_time).^2,1)));
        Tjv(:,:,1)       = bst_gain_orient(squeeze(Tjv(:,:,1))',GridOrient,GridAtlas);
        Tjv(:,:,2)       = bst_gain_orient(squeeze(Tjv(:,:,2))',GridOrient,GridAtlas);
        J_time_angle      = exp(i*angle(squeeze(Tjv(:,:,1))*V_time + squeeze(Tjv(:,:,2))*V_time));
        J_time           = J_time.*J_time_angle;
    end
    J_time               = J_time/(scaleKe(1)*scaleKe(2));
end
