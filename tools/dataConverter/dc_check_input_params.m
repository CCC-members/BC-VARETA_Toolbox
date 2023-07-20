function [status,errors,rejectedSub,dstruct] = dc_check_input_params(anatomyPath,functionalPath,outputPath,template)

status = true;
errors = [];
rejectedSub = [];
%% Checking Anatomy 
[~ ,~,ext] = fileparts(anatomyPath);
if(~isfolder(anatomyPath) && ~ismember(ext,{'.zip','.rar'}))
    errors = [errors; "The Anatomy path is not a folder or a zipped file"];
end
anatomies = dir(anatomyPath);
anatomies(ismember({anatomies.name},{'.','..'})) = [];
pos = 1;
for i=1:length(anatomies)
    anatomy = anatomies(i);
    subID = anatomy.name;
    if(~anatomy.isdir)
        subZip = fullfile(anatomy.folder,anatomy.name);
        [~,subID,~] = fileparts(subZip);
        unzip(subZip,fullfile(anatomy.folder,subID));
    end
    protocolFile = dir(fullfile(anatomy.folder,subID,'**','protocol.mat'));
    protocolInfo = load(fullfile(protocolFile.folder,protocolFile.name));
    subjects = protocolInfo.ProtocolSubjects.Subject; 
    studies = protocolInfo.ProtocolStudies.Study;
    for j=1:length(subjects)
        subject = subjects(j);
        subAnat = load(fullfile(anatomy.folder,subID,'anat',subject.FileName)); 
        dstruct(pos).Name = subject.Name;
        dstruct(pos).Anatomy.Head = fullfile(anatomy.folder,subID,'anat',subAnat.Scalp);
        dstruct(pos).Anatomy.Outerskull = fullfile(anatomy.folder,subID,'anat',subAnat.OuterSkull);
        dstruct(pos).Anatomy.Innerskull = fullfile(anatomy.folder,subID,'anat',subAnat.InnerSkull);
        dstruct(pos).Anatomy.Cortex = fullfile(anatomy.folder,subID,'anat',subAnat.Cortex);
        study = studies(contains({studies.FileName},subject.Name));
        study(ismember({study.Name},'@default_study')) = [];        
        dstruct(pos).Anatomy.Channel = fullfile(anatomy.folder,subID,'data',study.Channel.FileName);
        dstruct(pos).Modality = study.Channel.Modalities{1};
        dstruct(pos).Anatomy.HeadModel = fullfile(anatomy.folder,subID,'data',study.HeadModel.FileName);
        
        pos = pos + 1;
    end
end

%% Checking Functional 
if(~isfolder(functionalPath))
    errors = [errors; "The Functional path is not a folder"];
end
functionals = dir(functionalPath);
functionals(ismember({functionals.name},{'.','..'})) = [];
reject = [];
for i=1:length(dstruct)
    Name = dstruct(i).Name;
    fInfo = functionals(contains({functionals.name},Name));
    if (isempty(fInfo))
        reject = [reject, i];
    else
        dstruct(i).Functional = fullfile(fInfo.folder,fInfo.name);
    end
end
dstruct(reject) = [];

%% Checking Functional 
if(~isfolder(outputPath))
    errors = [errors; "The Output path is not a folder"];
    
end

