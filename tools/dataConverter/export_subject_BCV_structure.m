function export_error = export_subject_BCV_structure(properties,subID,CSurfaces,sub_to_FSAve)

%%
%% Get Protocol information
%%
export_error = [];

% Get subject directory
ProtocolInfo    = bst_get('ProtocolInfo');
sSubject        = bst_get('Subject', subID);
output_path     = properties.general_params.output_path;


% Get the current Study
[sStudies, ~]   = bst_get('StudyWithSubject', sSubject.FileName, 'intra_subject');
if(length(sStudies)>1)
    conditions  = [sStudies.Condition];
    sStudy      = sStudies(find(contains(conditions,strcat('@raw')),1));else
    sStudy      = sStudies;
end
if(isempty(sSubject) || isempty(sSubject.iAnatomy) || isempty(sSubject.iCortex) || isempty(sSubject.iInnerSkull) || isempty(sSubject.iOuterSkull) || isempty(sSubject.iScalp))
    return;
end

%% Uploding Subject file into BrainStorm Protocol
disp('BST-P ->> Uploading Subject file into BrainStorm Protocol.');
disp('==========================================================================');

%%
%% Genering leadfield file
%%
disp ("-->> Genering leadfield file");
[HeadModels,modality] = get_headmodels(ProtocolInfo.STUDIES,sStudy);

%%
%% Genering Channels file
%%
disp ("-->> Genering channels file");
if(isempty(sStudy.iChannel))
    sStudy.iChannel = 1;
end
BSTChannelsFile = bst_fullfile(ProtocolInfo.STUDIES,sStudy.Channel(sStudy.iChannel).FileName);
Cdata           = load(BSTChannelsFile);

%%
%% Getting surfaces
%%
[Shead, Sout, Sinn, Scortex] = get_surfaces(ProtocolInfo, sSubject, CSurfaces, sub_to_FSAve);

%%
%% Saving files in BC-VARETA Structure
%%
base_path = fullfile(output_path,'BC-V_Structure',ProtocolInfo.Comment);
save_output_files(base_path, modality, subID, HeadModels, Cdata, Shead, Sout, Sinn, Scortex);

end

