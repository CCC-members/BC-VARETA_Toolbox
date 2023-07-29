function  convertToBC_V_Structure(data,BC_V_Workspace)
for i=1:length(data)
    dataS                       = data(i);
    base_path                   = fullfile(BC_V_Workspace);
    modality                    = dataS.Modality;
    subID                       = dataS.Name;
    disp(strcat("-->> Exporting subject:" , subID));
    HeadModels.HeadModel        = load(dataS.Structural.HeadModel);
    Cdata                       = load(dataS.Structural.Channel);
    Shead                       = load(dataS.Structural.Head);
    Sout                        = load(dataS.Structural.Outerskull);
    Sinn                        = load(dataS.Structural.Innerskull);
    Scortex.Sc                  = load(dataS.Structural.Cortex);
    Scortex.iCortex             = 1;
    if(isequal(Scortex.Sc.Atlas(Scortex.Sc.iAtlas).Name,'Structures') || isempty(Scortex.Sc.Atlas(Scortex.Sc.iAtlas).Scouts))
        Scortex.Sc.Atlas.Name      = 'User scouts';
        Scortex.Sc.Atlas.Scouts    = generate_scouts(Scortex.Sc);
    end
    MEEG = load(dataS.Functional);
    action                      = 'new';
    if(isequal(modality,'EEG'))
        labels                  = {MEEG.chanlocs(:).labels};
    elseif(isequal(modality,'MEG'))
        labels                  = MEEG.labels;
    else
        labels                  = MEEG.dnames;
    end
    [Cdata, HeadModel]          = filter_structural_result_by_preproc_data(labels, Cdata, HeadModels.HeadModel);
    HeadModels.HeadModel        = HeadModel;
    HeadModels.iHeadModel       = 1;
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

