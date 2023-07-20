function [subject_report_path, report_name] = get_report_path(properties, subID)

%%
%% Checking the report output structure
%%
ProtocolName            = properties.general_params.bst_config.protocol_name;
report_output_path      = properties.general_params.reports.output_path;
if(report_output_path == "local")
    report_output_path = pwd;
end
if(~isfolder(report_output_path))
    mkdir(report_output_path);
end
if(~isfolder(fullfile(report_output_path,'Reports')))
    mkdir(fullfile(report_output_path,'Reports'));
end
if(~isfolder(fullfile(report_output_path,'Reports',ProtocolName)))
    mkdir(fullfile(report_output_path,'Reports',ProtocolName));
end
if(~isfolder(fullfile(report_output_path,'Reports',ProtocolName,subID)))
    mkdir(fullfile(report_output_path,'Reports',ProtocolName,subID));
end
subject_report_path = fullfile(report_output_path,'Reports',ProtocolName,subID);
report_name         = fullfile(subject_report_path,[subID,'.html']);
iter                = 2;
while(isfile(report_name))
    report_name     = fullfile(subject_report_path,[subID,'_Iter_', num2str(iter),'.html']);
    iter            = iter + 1;
end

end

