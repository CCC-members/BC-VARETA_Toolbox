function  convertToBC_V_Structure(data,BC_V_Workspace)
for i=1:length(data)
    dataS = data(i);
    base_path = fullfile(BC_V_Workspace);
    modality = dataS.Modality;
    subID = dataS.Name;
    disp(strcat("-->> Exporting subject:" , subID));
    HeadModels.HeadModel = load(dataS.Structural.HeadModel);
    Cdata = load(dataS.Structural.Channel);
    Shead = load(dataS.Structural.Head);
    Sout = load(dataS.Structural.Outerskull);
    Sinn = load(dataS.Structural.Innerskull);
    Scortex = load(dataS.Structural.Cortex);
    if(isequal(Scortex.Atlas(Scortex.iAtlas).Name,'Structures') || isempty(Scortex.Atlas(Scortex.iAtlas).Scouts))
        Scortex.Atlas.Name    = 'User scouts';
        Scortex.Atlas.Scouts  = generate_scouts(Cortex);
    end
    MEEG = load(dataS.Functional);
    action = 'new';
    if(isequal(modality,'EEG'))
        labels              = {MEEG.chanlocs(:).labels};
    elseif(isequal(modality,'MEG'))
        labels              = MEEG.labels;
    else
        labels              = MEEG.dnames;
    end
    [Cdata, HeadModel]      = filter_structural_result_by_preproc_data(labels, Cdata, HeadModels.HeadModel);
    HeadModels.HeadModel    = HeadModel;
    save_output_files(action, ...
        base_path, ...
        modality, ...
        subID, ...
        HeadModels, ...
        Cdata, ...
        Shead, ...
        Sout, ...
        Sinn, ...
        Scortex, ...
        MEEG);
end
end

