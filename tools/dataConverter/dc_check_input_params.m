function [status,errors,rejectedSub,data] = dc_check_input_params(structuralPath,functionalPath,outputPath)

status = true;
errors = [];
rejectedSub = {};
data = [];

%% Checking Functional
if(~isfolder(outputPath))
    errors = "The Output path is not a folder";
    status = false;
    return;
end

%% Checking Anatomy
[~ ,~,ext] = fileparts(structuralPath);
if(~isfolder(structuralPath) && ~ismember(ext,{'.zip','.rar'}))
    errors = [errors; "The Structural path is not a folder or a zipped file"];
    status = false;
    return;
end

%% Checking Functional
if(~isfolder(functionalPath))
    errors = [errors; "The Functional path is not a folder"];
end

%% Defining cases
folders = dir(structuralPath);
folders(ismember({folders.name},{'.','..'})) = [];
if(ismember({folders.name},{'anat','data'}))
    structure = 'Protocol';
else
    structure = 'Subjects';
end
% structure = 'Template';

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
                data(pos).Name = subject.Name;
                data(pos).Structural.Head = fullfile(structuralPath,'anat',subAnat.Scalp);
                data(pos).Structural.Outerskull = fullfile(structuralPath,'anat',subAnat.OuterSkull);
                data(pos).Structural.Innerskull = fullfile(structuralPath,'anat',subAnat.InnerSkull);
                data(pos).Structural.Cortex = fullfile(structuralPath,'anat',subAnat.Cortex);
                study = studies(contains({studies.FileName},subject.Name));
                study(ismember({study.Name},'@default_study')) = [];
                data(pos).Structural.Channel = fullfile(structuralPath,'data',study.Channel.FileName);
                data(pos).Modality = study.Channel.Modalities{1};
                data(pos).Structural.HeadModel = fullfile(structuralPath,'data',study.HeadModel.FileName);
                pos = pos + 1;
            catch ME
                rejectedSub = [rejectedSub; subject.Name];
                errors = [errors; ME.message];
            end
        end
        functionals = dir(functionalPath);
        functionals(ismember({functionals.name},{'.','..'})) = [];
        reject = [];        
        for i=1:length(data)
            Name = data(i).Name;
            fInfo = functionals(contains({functionals.name},Name));
            if (isempty(fInfo))
                reject = [reject, i];
                errors = [errors; "Functional data not found"];
            else
                data(i).Functional = fullfile(fInfo.folder,fInfo.name);
            end
        end
        rejectedSub = data(reject);
        data(reject) = [];
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
            reject = [];
            for j=1:length(subjects)
                try
                    subject = subjects(j);
                    subAnat = load(fullfile(anatomy.folder,subID,'anat',subject.FileName));
                    data(pos).Name = subject.Name;
                    data(pos).Structural.Head = fullfile(anatomy.folder,subID,'anat',subAnat.Scalp);
                    data(pos).Structural.Outerskull = fullfile(anatomy.folder,subID,'anat',subAnat.OuterSkull);
                    data(pos).Structural.Innerskull = fullfile(anatomy.folder,subID,'anat',subAnat.InnerSkull);
                    data(pos).Structural.Cortex = fullfile(anatomy.folder,subID,'anat',subAnat.Cortex);
                    study = studies(contains({studies.FileName},subject.Name));
                    study(ismember({study.Name},'@default_study')) = [];
                    data(pos).Structural.Channel = fullfile(anatomy.folder,subID,'data',study.Channel.FileName);
                    data(pos).Modality = study.Channel.Modalities{1};
                    data(pos).Structural.HeadModel = fullfile(anatomy.folder,subID,'data',study.HeadModel.FileName);
                    pos = pos + 1;
                catch ME
                    reject = [reject, i];
                    errors = [errors; ME.message];
                end
            end
        end
        rejectedSub = subjects(reject);
        functionals = dir(functionalPath);
        functionals(ismember({functionals.name},{'.','..'})) = [];
        reject = [];
        for i=1:length(data)
            Name = data(i).Name;
            fInfo = functionals(contains({functionals.name},Name));
            if (isempty(fInfo))
                reject = [reject, i];
                errors = [errors; "Functional data not found"];
            else
                data(i).Functional = fullfile(fInfo.folder,fInfo.name);
            end
        end
        rejectedSub = subjects(reject);
        data(reject) = [];
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
            [~,data(pos).Name,~] = fileparts(functionals(i).name);
            data(pos).Functional = fInfo;
            data(pos).Structural.Head = fullfile(anatomy.folder,templateID,'anat',subAnat.Scalp);
            data(pos).Structural.Outerskull = fullfile(anatomy.folder,templateID,'anat',subAnat.OuterSkull);
            data(pos).Structural.Innerskull = fullfile(anatomy.folder,templateID,'anat',subAnat.InnerSkull);
            data(pos).Structural.Cortex = fullfile(anatomy.folder,templateID,'anat',subAnat.Cortex);
            data(pos).Structural.Channel = fullfile(anatomy.folder,templateID,'data',study.Channel.FileName);
            data(pos).Modality = study.Channel.Modalities{1};
            data(pos).Structural.HeadModel = fullfile(anatomy.folder,templateID,'data',study.HeadModel.FileName);
            pos = pos + 1;
        end
end
end

