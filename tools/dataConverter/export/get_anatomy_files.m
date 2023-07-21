function [outputArg1,outputArg2] = get(anatomy,inputArg2)

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

