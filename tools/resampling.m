clear all
%% data directoris containing the files to interpolate
data_dir11 = 'E:\BC-V_Activation1\FSAverage\MEG\eLORETA';
data_dir12 = 'E:\BC-V_Activation1\FSAverage\MEG\LCMV';
data_dir13 = 'E:\BC-V_Activation1\FSAverage\MEG\sSSBLpp';
data_dir21 = 'E:\BC-V_Activation1\Native\MEG\eLORETA';
data_dir22 = 'E:\BC-V_Activation1\Native\MEG\LCMV';
data_dir23 = 'E:\BC-V_Activation1\Native\MEG\sSSBLpp';
data_dir   = {data_dir11,data_dir12,data_dir13,data_dir21,data_dir22,data_dir23};

%% interpolation space
ratio   = 508/408; % frequency dilatation ratio
fspace  = 0.1:0.5:(0.1+99*0.5); % frequency space to interpolate
fspace0 = fspace*ratio; % actual frequency space of the data

%% files to interpolate
var_name = {'F3Dstat','J3D','J3Dnorm','J3Dnorm_sp','J3Dsp'};
for i = 1:length(data_dir)
    current_dir = data_dir{i};
    for ii = 1:5
        file_name{ii}  = strcat(current_dir,'\',var_name{ii});
        load(file_name{ii});
    end
    %% cycle for sources and subjects
    for subject = 1:70
        for source = 1:8002
            % F3Dstat
            act_vector                   = squeeze(F3Dstat(source,:,subject));
            act_vector_interpolated      = interp1(fspace0,act_vector,fspace,'spline');
            F3Dstat(source,:,subject)    = act_vector_interpolated;
            % J3D
            act_vector                   = squeeze(J3D(source,:,subject));
            act_vector_interpolated      = interp1(fspace0,act_vector,fspace,'spline');
            J3D(source,:,subject)        = act_vector_interpolated;
            % J3Dnorm
            act_vector                   = squeeze(J3Dnorm(source,:,subject));
            act_vector_interpolated      = interp1(fspace0,act_vector,fspace,'spline');
            J3Dnorm(source,:,subject)    = act_vector_interpolated;
            % J3Dnorm_sp
            act_vector                   = squeeze(J3Dnorm_sp(source,:,subject));
            act_vector_interpolated      = interp1(fspace0,act_vector,fspace,'spline');
            J3Dnorm_sp(source,:,subject) = act_vector_interpolated;
            % J3Dsp
            act_vector                   = squeeze(J3Dsp(source,:,subject));
            act_vector_interpolated      = interp1(fspace0,act_vector,fspace,'spline');
            J3Dsp(source,:,subject)      = act_vector_interpolated;
        end
    end
    %% saving variables
    for ii = 1:5
        interpolated_var_name = strcat(file_name{ii},'_interp');
        save(interpolated_var_name,var_name{ii});
    end
end
