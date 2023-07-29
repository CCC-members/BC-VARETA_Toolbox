function save_error = save_output_files(varargin)

save_error = [];
%%
%% Creating structure 
%%

for i=1:length(varargin)
   eval([inputname(i) '= varargin{i};']); 
end
if(isequal(action,'new'))
    disp(strcat("-->> Creating subject output structure"));
    action = 'all';
    [output_subject_dir]        = create_data_structure(base_path,subID,action);   
    subject_info                = struct;
    subject_info.name           = subID;
    subject_info.modality       = modality;
    dirref                      = replace(fullfile('meeg','meeg.mat'),'\','/');
    subject_info.meeg_dir       = dirref;
    dirref                      = replace(fullfile('leadfield','headmodel.mat'),'\','/');
    subject_info.leadfield_dir  = dirref;
    dirref                      = replace(fullfile('surf','surf.mat'),'\','/');
    subject_info.surf_dir       = dirref;
    dirref                      = replace(fullfile('channel','channel.mat'),'\','/');
    subject_info.channel_dir    = dirref;
    dirref                      = replace(fullfile('scalp','scalp.mat'),'\','/');
    subject_info.scalp_dir      = dirref;
    dirref                      = replace(fullfile('scalp','innerskull.mat'),'\','/');
    subject_info.innerskull_dir = dirref;
    dirref                      = replace(fullfile('scalp','outerskull.mat'),'\','/');
    subject_info.outerskull_dir = dirref;
    subject_info.completed      = true;
    
    % Saving subject files
    disp ("-->> Saving MEEG file");
    save(fullfile(output_subject_dir,subject_info.meeg_dir),'-struct','MEEG');
    disp ("-->> Saving channel file");
    save(fullfile(output_subject_dir,subject_info.channel_dir),'-struct','Cdata');
    disp ("-->> Saving leadfield file");
    save(fullfile(output_subject_dir,subject_info.leadfield_dir),'-struct','HeadModels');
    disp ("-->> Saving scalp file");
    save(fullfile(output_subject_dir,'scalp','scalp.mat'),'-struct','Shead');
    disp ("-->> Saving outer skull file");
    save(fullfile(output_subject_dir,'scalp','outerskull.mat'),'-struct','Sout');
    disp ("-->> Saving inner skull file");
    save(fullfile(output_subject_dir,'scalp','innerskull.mat'),'-struct','Sinn');
    disp ("-->> Saving surf file");
    save(fullfile(output_subject_dir,'surf','surf.mat'),'-struct','Scortex');
    disp ("-->> Saving subject file");
    save(fullfile(output_subject_dir,'subject.mat'),'-struct','subject_info');
end
if(isequal(action,'update'))
    % Updating subject files
    dirref                  = replace(fullfile('meeg','meeg.mat'),'\','/');
    subject_info.meeg_dir   = dirref;
    subject_info.completed  = true;
    
    disp ("-->> Saving MEEG file");
    if(~isfolder(fullfile(base_path,subID,'meeg')))
        mkdir(fullfile(base_path,subID,'meeg'));
    end
    save(fullfile(base_path,subID,subject_info.meeg_dir),'-struct','MEEG');
    disp ("-->> Saving channel file");
    save(fullfile(base_path,subID,subject_info.channel_dir),'-struct','Cdata');
    disp ("-->> Saving leadfield file");
    save(fullfile(base_path,subID,subject_info.leadfield_dir),'-struct','HeadModels');
    disp ("-->> Saving subject file");
    save(fullfile(base_path,subID,'subject.mat'),'-struct','subject_info');
end
if(isequal(action,'event'))
    % Creating subject folder structure
    disp(strcat("-->> Creating subject output structure"));
    base_path = fullfile(base_path,MEEG.event_name);
    if(~isfolder(base_path))
        mkdir(base_path);
    end
    action = 'all';
    [output_subject_dir]        = create_data_structure(base_path,subID,action);        
    subject_info.name           = subID;
    dirref                      = replace(fullfile('meeg','meeg.mat'),'\','/');
    subject_info.meeg_dir       = dirref;
    dirref                      = replace(fullfile('leadfield','headmodel.mat'),'\','/');
    subject_info.leadfield_dir  = dirref;
    dirref                      = replace(fullfile('surf','surf.mat'),'\','/');
    subject_info.surf_dir       = dirref;
    dirref                      = replace(fullfile('channel','channel.mat'),'\','/');
    subject_info.channel_dir    = dirref;
    dirref                      = replace(fullfile('scalp','scalp.mat'),'\','/');
    subject_info.scalp_dir      = dirref;
    dirref                      = replace(fullfile('scalp','innerskull.mat'),'\','/');
    subject_info.innerskull_dir = dirref;
    dirref                      = replace(fullfile('scalp','outerskull.mat'),'\','/');
    subject_info.outerskull_dir = dirref;
    subject_info.event          = MEEG.event_name;
    subject_info.completed      = true;
    
    % Saving subject files
    disp ("-->> Saving MEEG file");
    save(fullfile(output_subject_dir,subject_info.meeg_dir),'-struct','MEEG');
    disp ("-->> Saving channel file");
    save(fullfile(output_subject_dir,subject_info.channel_dir),'-struct','Cdata');
    disp ("-->> Saving leadfield file");
    save(fullfile(output_subject_dir,subject_info.leadfield_dir),'-struct','HeadModels');
    disp ("-->> Saving scalp file");
    save(fullfile(output_subject_dir,'scalp','scalp.mat'),'-struct','Shead');
    disp ("-->> Saving outer skull file");
    save(fullfile(output_subject_dir,'scalp','outerskull.mat'),'-struct','Sout');
    disp ("-->> Saving inner skull file");
    save(fullfile(output_subject_dir,'scalp','innerskull.mat'),'-struct','Sinn');
    disp ("-->> Saving surf file");
    save(fullfile(output_subject_dir,'surf','surf.mat'),'-struct','Scortex');
    disp ("-->> Saving subject file");
    save(fullfile(output_subject_dir,'subject.mat'),'-struct','subject_info');    
end


end



