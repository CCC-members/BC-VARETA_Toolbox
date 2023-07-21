function [status,errors,rejectedSub,dstruct] = dc_check_input_params(structuralPath,functionalPath,outputPath,varargin)

for k = 4 : nargin
    eval([inputname(k) '=  varargin{k-3};']);
end

status = true;
errors = [];
rejectedSub = [];

%% Checking Functional
if(~isfolder(outputPath))
    errors = "The Output path is not a folder";
    status = false;
    return;
end

%% Checking Anatomy
[~ ,~,ext] = fileparts(structuralPath);
if(~isfolder(structuralPath) && ~ismember(ext,{'.zip','.rar'}))
    errors = [errors; "The Anatomy path is not a folder or a zipped file"];
    status = false;
    return;
end

%% Checking Functional
if(~isfolder(functionalPath))
    errors = [errors; "The Functional path is not a folder"];
end

switch structure
    case 'Protocol'
        protocolInfo = load(fullfile(structuralPath,'data','protocol.mat'));
        subjects = protocolInfo.ProtocolSubjects.Subject;
        studies = protocolInfo.ProtocolStudies.Study;
        pos = 1;
        for j=1:length(subjects)
            subject = subjects(j);
            try                
                subAnat = load(fullfile(structuralPath,'anat',subject.FileName));
                dstruct(pos).Name = subject.Name;
                dstruct(pos).Anatomy.Head = fullfile(structuralPath,'anat',subAnat.Scalp);
                dstruct(pos).Anatomy.Outerskull = fullfile(structuralPath,'anat',subAnat.OuterSkull);
                dstruct(pos).Anatomy.Innerskull = fullfile(structuralPath,'anat',subAnat.InnerSkull);
                dstruct(pos).Anatomy.Cortex = fullfile(structuralPath,'anat',subAnat.Cortex);
                study = studies(contains({studies.FileName},subject.Name));
                study(ismember({study.Name},'@default_study')) = [];
                dstruct(pos).Anatomy.Channel = fullfile(structuralPath,'data',study.Channel.FileName);
                dstruct(pos).Modality = study.Channel.Modalities{1};
                dstruct(pos).Anatomy.HeadModel = fullfile(structuralPath,'data',study.HeadModel.FileName);
                pos = pos + 1;
            catch ME
                rejectedSub = [rejectedSub; subject.Name];
                errors = [errors; ME.message];
            end
        end
        functionals = dir(functionalPath);
        functionals(ismember({functionals.name},{'.','..'})) = [];
        reject = [];
        for i=1:length(dstruct)
            Name = dstruct(i).Name;
            fInfo = functionals(contains({functionals.name},Name));
            if (isempty(fInfo))
                reject = [reject, i];
                rejectedSub = [rejectedSub; Name];
                errors = [errors; "Functional data not found"];
            else
                dstruct(i).Functional = fullfile(fInfo.folder,fInfo.name);
            end
        end
        dstruct(reject) = [];
    case 'Subjects'
        anatomies = dir(structuralPath);
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
                try
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
                catch ME
                    rejectedSub = [rejectedSub; subject.Name];
                    errors = [errors; ME.message];
                end
            end
        end
        functionals = dir(functionalPath);
        functionals(ismember({functionals.name},{'.','..'})) = [];
        reject = [];
        for i=1:length(dstruct)
            Name = dstruct(i).Name;
            fInfo = functionals(contains({functionals.name},Name));
            if (isempty(fInfo))
                reject = [reject, i];
                rejectedSub = [rejectedSub; Name];
                errors = [errors; "Functional data not found"];
            else
                dstruct(i).Functional = fullfile(fInfo.folder,fInfo.name);
            end
        end
        dstruct(reject) = [];
    case 'Template'
        anatomies = dir(structuralPath);
        anatomies(ismember({anatomies.name},{'.','..'})) = [];
        anatomies(~contains({anatomies.name},templateName)) = [];
        if(isempty(anatomies))
            status = false;
            errors = "Template name does not found";
            return;
        end
        anatomy = anatomies;
        if(~anatomy.isdir)
            subZip = fullfile(anatomy.folder,anatomy.name);
            [~,templateID,~] = fileparts(subZip);
            unzip(subZip,fullfile(anatomy.folder,templateID));
        end
        protocolFile = dir(fullfile(anatomy.folder,templateID,'**','protocol.mat'));
        protocolInfo = load(fullfile(protocolFile.folder,protocolFile.name));
        subject = protocolInfo.ProtocolSubjects.Subject;
        subAnat = load(fullfile(anatomy.folder,templateID,'anat',subject.FileName));

        studies = protocolInfo.ProtocolStudies.Study;
        study = studies(contains({studies.FileName},subject.Name));
        study(ismember({study.Name},'@default_study')) = [];
        functionals = dir(functionalPath);
        functionals(ismember({functionals.name},{'.','..'})) = [];
        pos = 1;
        for i=1:length(functionals)
            fInfo = fullfile(functionals(i).folder,functionals(i).name);
            [~,dstruct(pos).Name,~] = fileparts(functionals(i).name);
            dstruct(pos).Functional = fInfo;
            dstruct(pos).Anatomy.Head = fullfile(anatomy.folder,templateID,'anat',subAnat.Scalp);
            dstruct(pos).Anatomy.Outerskull = fullfile(anatomy.folder,templateID,'anat',subAnat.OuterSkull);
            dstruct(pos).Anatomy.Innerskull = fullfile(anatomy.folder,templateID,'anat',subAnat.InnerSkull);
            dstruct(pos).Anatomy.Cortex = fullfile(anatomy.folder,templateID,'anat',subAnat.Cortex);
            dstruct(pos).Anatomy.Channel = fullfile(anatomy.folder,templateID,'data',study.Channel.FileName);
            dstruct(pos).Modality = study.Channel.Modalities{1};
            dstruct(pos).Anatomy.HeadModel = fullfile(anatomy.folder,templateID,'data',study.HeadModel.FileName);
            pos = pos + 1;
        end
end
end

