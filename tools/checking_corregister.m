path = 'E:\Data\CHBM\HeadModel_CHBM_0\anat\';
subjects = dir(path);
compare = struct;
for j=4:length(subjects)
    subject = subjects(j);
    subject_name = subject.name;
    disp(strcat("-->> ", subject.name));
    cortex_8k = load(fullfile(subject.folder,subject.name,'tess_cortex_concat.mat'));
    cortex_32k = load(fullfile(subject.folder,subject.name,'tess_cortex_concat_02.mat'));
    ind = zeros(length(cortex_8k.Vertices),1);
    x1  = cortex_32k.Vertices;
    for i = 1:length(cortex_8k.Vertices)
        x2     = repmat(cortex_8k.Vertices(i,:),length(x1),1);
        dist   = sqrt(sum(abs(x2-x1).^2,2));
        ind(i) = find(dist == 0);
    end
    compare(j).name = subject_name;
    compare(j).ind  = ind;
end
            

