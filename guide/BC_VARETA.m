classdef BC_VARETA < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        BCVARETAToolboxv10UIFigure     matlab.ui.Figure
        FileMenu                       matlab.ui.container.Menu
        NewprocessMenu                 matlab.ui.container.Menu
        TestdataMenu                   matlab.ui.container.Menu
        TestDataDownloadMenu           matlab.ui.container.Menu
        TestDataImportMenu             matlab.ui.container.Menu
        TestDataRunMenu                matlab.ui.container.Menu
        DatasetMenu                    matlab.ui.container.Menu
        DatasetLoadMenu                matlab.ui.container.Menu
        DatasetListMenu                matlab.ui.container.Menu
        DatasetDeleteMenu              matlab.ui.container.Menu
        ExitMenu                       matlab.ui.container.Menu
        ToolsMenu                      matlab.ui.container.Menu
        GroupstastMenu                 matlab.ui.container.Menu
        ViewMenu                       matlab.ui.container.Menu
        DatasetresultsMenu             matlab.ui.container.Menu
        HelpMenu                       matlab.ui.container.Menu
        CheckforupdateMenu             matlab.ui.container.Menu
        DocumentationMenu              matlab.ui.container.Menu
        AboutMenu                      matlab.ui.container.Menu
        Toolbar                        matlab.ui.container.Toolbar
        NewProcessPushTool             matlab.ui.container.toolbar.PushTool
        Default_ParamsPushTool         matlab.ui.container.toolbar.PushTool
        ListDatasetPushTool            matlab.ui.container.toolbar.PushTool
        CancelButton                   matlab.ui.control.Button
        RunButton                      matlab.ui.control.Button
        TabGroup                       matlab.ui.container.TabGroup
        PrincipalTab                   matlab.ui.container.Tab
        GeneralTab                     matlab.ui.container.Tab
        gener_PrincipalParamsPanel     matlab.ui.container.Panel
        gener_DescriptionTextArea      matlab.ui.control.TextArea
        DescriptionTextAreaLabel       matlab.ui.control.Label
        InputfolderLabel_2             matlab.ui.control.Label
        gener_DatasetNameEditField     matlab.ui.control.EditField
        gener_RuninGPUSwitch           matlab.ui.control.Switch
        RuninGPUSwitchLabel            matlab.ui.control.Label
        GetsystemresponseSwitchLabel   matlab.ui.control.Label
        gener_GetsystemresponseSwitch  matlab.ui.control.Switch
        RunbyfrequencybinSwitchLabel   matlab.ui.control.Label
        gener_RunbyfrequencybinSwitch  matlab.ui.control.Switch
        RunbytrialsSwitchLabel         matlab.ui.control.Label
        gener_RunbytrialsSwitch        matlab.ui.control.Switch
        AnalysislevelDropDownLabel     matlab.ui.control.Label
        gener_AnalysislevelDropDown    matlab.ui.control.DropDown
        OutputfolderLabel              matlab.ui.control.Label
        gener_OutputfolderEditField    matlab.ui.control.EditField
        InputfolderLabel               matlab.ui.control.Label
        gener_InputfolderEditField     matlab.ui.control.EditField
        gener_SelectInputButton        matlab.ui.control.Button
        gener_SelectOutputButton       matlab.ui.control.Button
        SensorTab                      matlab.ui.container.Tab
        spect_OtherparamsPanel         matlab.ui.container.Panel
        sens_FFTMethodDropDown         matlab.ui.control.DropDown
        FTTmethodLabel                 matlab.ui.control.Label
        spect_win_orderSlider          matlab.ui.control.Slider
        win_orderLabel                 matlab.ui.control.Label
        spect_freq_gfiltvarEditField   matlab.ui.control.NumericEditField
        freq_gfiltvarEditFieldLabel    matlab.ui.control.Label
        spect_FreqmaxEditField         matlab.ui.control.NumericEditField
        FreqmaxEditFieldLabel          matlab.ui.control.Label
        spect_sample_freqEditField     matlab.ui.control.NumericEditField
        SamplefreqLabel                matlab.ui.control.Label
        spect_freq_resolEditField      matlab.ui.control.NumericEditField
        FreqresolutionLabel            matlab.ui.control.Label
        spect_FreqPanel                matlab.ui.container.Panel
        GridLayout                     matlab.ui.container.GridLayout
        Gamma2Label                    matlab.ui.control.Label
        Gamma1Label                    matlab.ui.control.Label
        BetaLabel                      matlab.ui.control.Label
        AlphaLabel                     matlab.ui.control.Label
        ThetaLabel                     matlab.ui.control.Label
        DeltaLabel                     matlab.ui.control.Label
        gamma2_runSwitch               matlab.ui.control.Switch
        gamma1_runSwitch               matlab.ui.control.Switch
        beta_runSwitch                 matlab.ui.control.Switch
        alpha_runSwitch                matlab.ui.control.Switch
        theta_runSwitch                matlab.ui.control.Switch
        delta_runSwitch                matlab.ui.control.Switch
        gamma2_endEditField            matlab.ui.control.NumericEditField
        gamma2_startEditField          matlab.ui.control.NumericEditField
        gamma1_endEditField            matlab.ui.control.NumericEditField
        gamma1_startEditField          matlab.ui.control.NumericEditField
        beta_endEditField              matlab.ui.control.NumericEditField
        beta_startEditField            matlab.ui.control.NumericEditField
        alpha_endEditField             matlab.ui.control.NumericEditField
        alpha_startEditField           matlab.ui.control.NumericEditField
        theta_endEditField             matlab.ui.control.NumericEditField
        theta_startEditField           matlab.ui.control.NumericEditField
        delta_endEditField             matlab.ui.control.NumericEditField
        delta_startEditField           matlab.ui.control.NumericEditField
        EndLabel                       matlab.ui.control.Label
        StartLabel                     matlab.ui.control.Label
        ProcessLabel                   matlab.ui.control.Label
        NameLabel                      matlab.ui.control.Label
        ActivationTab                  matlab.ui.container.Tab
        activ_OtherparamsPanel         matlab.ui.container.Panel
        FieldLabel                     matlab.ui.control.Label
        activ_fieldSlider              matlab.ui.control.Slider
        regLaplacianEditFieldLabel     matlab.ui.control.Label
        activ_regLaplacianEditField    matlab.ui.control.NumericEditField
        bSulcEditFieldLabel            matlab.ui.control.Label
        activ_bSulcEditField           matlab.ui.control.NumericEditField
        bGiriEditFieldLabel            matlab.ui.control.Label
        activ_bGiriEditField           matlab.ui.control.NumericEditField
        aSulcEditFieldLabel            matlab.ui.control.Label
        activ_aSulcEditField           matlab.ui.control.NumericEditField
        aGiriLabel                     matlab.ui.control.Label
        activ_aGiriEditField           matlab.ui.control.NumericEditField
        activ_IsCurvCheckBox           matlab.ui.control.CheckBox
        activ_IsParcelCheckBox         matlab.ui.control.CheckBox
        activ_IsNeighCheckBox          matlab.ui.control.CheckBox
        IsCurvLabel                    matlab.ui.control.Label
        IsParcelLabel                  matlab.ui.control.Label
        IsNeighLabel                   matlab.ui.control.Label
        activ_MethodsPanel             matlab.ui.container.Panel
        lcmv_deltagammaEditField       matlab.ui.control.NumericEditField
        deltagammaEditField_2Label     matlab.ui.control.Label
        lcmv_gamma2EditField           matlab.ui.control.NumericEditField
        gamma2EditField_2Label         matlab.ui.control.Label
        lcmv_gamma1EditField           matlab.ui.control.NumericEditField
        gamma1EditField_2Label         matlab.ui.control.Label
        lcmv_thEditField               matlab.ui.control.NumericEditField
        ThresholdEditField_3Label      matlab.ui.control.Label
        eloreta_deltagammaEditField    matlab.ui.control.NumericEditField
        deltagammaEditFieldLabel       matlab.ui.control.Label
        eloreta_gamma2EditField        matlab.ui.control.NumericEditField
        gamma2EditFieldLabel           matlab.ui.control.Label
        eloreta_gamma1EditField        matlab.ui.control.NumericEditField
        gamma1EditFieldLabel           matlab.ui.control.Label
        eloreta_thEditField            matlab.ui.control.NumericEditField
        ThresholdEditField_2Label      matlab.ui.control.Label
        sssblpp_thEditField            matlab.ui.control.NumericEditField
        ThresholdEditFieldLabel        matlab.ui.control.Label
        CheckBox_lcmv                  matlab.ui.control.CheckBox
        LCMVLabel                      matlab.ui.control.Label
        CheckBox_eloreta               matlab.ui.control.CheckBox
        eLORETALabel                   matlab.ui.control.Label
        CheckBox_sssblpp               matlab.ui.control.CheckBox
        SSSBLppLabel                   matlab.ui.control.Label
        ConnectivityTab                matlab.ui.container.Tab
        conn_OtherparamsPanel          matlab.ui.container.Panel
        prewarmingLabel                matlab.ui.control.Label
        IsNeighLabel_2                 matlab.ui.control.Label
        IsCurvLabel_2                  matlab.ui.control.Label
        conn_fieldSlider               matlab.ui.control.Slider
        FieldLabel_2                   matlab.ui.control.Label
        conn_eigregEditField           matlab.ui.control.NumericEditField
        eigregLabel                    matlab.ui.control.Label
        conn_rth2EditField             matlab.ui.control.NumericEditField
        rth2Label                      matlab.ui.control.Label
        conn_rth1EditField             matlab.ui.control.NumericEditField
        rth1Label                      matlab.ui.control.Label
        conn_penaltySlider             matlab.ui.control.Slider
        penaltySliderLabel             matlab.ui.control.Label
        conn_prewarmingCheckBox        matlab.ui.control.CheckBox
        conn_axiEditField              matlab.ui.control.NumericEditField
        axiLabel                       matlab.ui.control.Label
        conn_ntrySlider                matlab.ui.control.Slider
        ntryLabel                      matlab.ui.control.Label
        conn_MaxiterinnerEditField     matlab.ui.control.NumericEditField
        MaxiterinnerEditFieldLabel     matlab.ui.control.Label
        conn_MaxiterouterEditField     matlab.ui.control.NumericEditField
        MaxiterouterEditFieldLabel     matlab.ui.control.Label
        conn_regLaplacianEditField     matlab.ui.control.NumericEditField
        regLaplacianEditField_2Label   matlab.ui.control.Label
        conn_bSulcEditField            matlab.ui.control.NumericEditField
        bSulcEditField_2Label          matlab.ui.control.Label
        conn_bGiriEditField            matlab.ui.control.NumericEditField
        bGiriEditField_2Label          matlab.ui.control.Label
        conn_aSulcEditField            matlab.ui.control.NumericEditField
        aSulcEditField_2Label          matlab.ui.control.Label
        conn_aGiriEditField            matlab.ui.control.NumericEditField
        aGiriEditField_2Label          matlab.ui.control.Label
        conn_IsNeighCheckBox           matlab.ui.control.CheckBox
        conn_IsCurvCheckBox            matlab.ui.control.CheckBox
        conn_MethodsPanel              matlab.ui.container.Panel
        hg_lasso_thEditField           matlab.ui.control.NumericEditField
        ThresholdEditField_2Label_2    matlab.ui.control.Label
        CheckBox_hg_lasso              matlab.ui.control.CheckBox
        hg_LASSOLabel                  matlab.ui.control.Label
        higgs_thEditField              matlab.ui.control.NumericEditField
        ThresholdEditFieldLabel_2      matlab.ui.control.Label
        CheckBox_higgs                 matlab.ui.control.CheckBox
        HIGGSLabel                     matlab.ui.control.Label
    end

    
    properties (Access = private)
        bcv_properties
        app_params; % Description        
        Name;
        BCVdir;
        Datasets_file;
        Datasets;
    end
    
    properties (Access = public)        
         % Description
    end
    
    methods (Access = private)
        
        function setPromptFcn(app,jTextArea)
            % Prevent overlapping reentry due to prompt replacement
            persistent inProgress
            if isempty(inProgress)
                inProgress = 1;  %#ok unused
            else
                return;
            end
            
            try
                % *** Prompt modification code goes here ***
                cwText = char(jTextArea.getText);
                app.gener_DescriptionTextArea.Value = cwText;
                % force prompt-change callback to fizzle-out...
                pause(0.02);
            catch
                % Never mind - ignore errors...
            end
            
            % Enable new callbacks now that the prompt has been modified
            inProgress = [];
            
        end  % setPromptFcn
        
        % Getting default properties and seting on user interface
        function set_default_properties(app)         
            app.bcv_properties.general_params          = app.bcv_properties.defaults.general_params;
            app.bcv_properties.sensor_params           = app.bcv_properties.defaults.sensor_params;
            app.bcv_properties.activation_params       = app.bcv_properties.defaults.activation_params;
            app.bcv_properties.connectivity_params     = app.bcv_properties.defaults.connectivity_params;
            % Saving property files
            saveJSON(app.bcv_properties.general_params,'bcv_properties/general_params.json');
            h = matlab.desktop.editor.openDocument(fullfile(pwd,'bcv_properties/general_params.json'));
            h.smartIndentContents
            h.save
            h.close
            saveJSON(app.bcv_properties.sensor_params,'bcv_properties/sensor_params.json');
            h = matlab.desktop.editor.openDocument(fullfile(pwd,'bcv_properties/sensor_params.json'));
            h.smartIndentContents
            h.save
            h.close
            saveJSON(app.bcv_properties.activation_params,'bcv_properties/activation_params.json');
            h = matlab.desktop.editor.openDocument(fullfile(pwd,'bcv_properties/activation_params.json'));
            h.smartIndentContents
            h.save
            h.close
            saveJSON(app.bcv_properties.connectivity_params,'bcv_properties/connectivity_params.json');
            h = matlab.desktop.editor.openDocument(fullfile(pwd,'bcv_properties/connectivity_params.json'));
            h.smartIndentContents
            h.save
            h.close
            
            %% Set the component values from properties files
            set_component_values(app);

        end % set_default_properties
        
        % Checking user interface params
        function checked = check_user_params(app)
            checked = true;
            if(isempty(app.gener_DatasetNameEditField.Value))
                msgbox('The Dataset name can not be empty.','Info');
                checked = false;
            end
            if(isempty(app.gener_InputfolderEditField.Value) || isempty(app.gener_OutputfolderEditField.Value))
                msgbox('The Input folder and Output folder fields can not be empty.','Info');
                checked = false;
            end
            if(~isfolder(app.gener_InputfolderEditField.Value))
                msgbox('The Input path is not a real folder.','Info');
                checked = false;
            end
            if(~isfolder(app.gener_OutputfolderEditField.Value))
                msgbox('The Output path is not a real folder.','Info');
                checked = false;
            else
                [~,values] = fileattrib(app.gener_OutputfolderEditField.Value);
                if(~values.UserWrite)
                    fprintf(2,strcat('\nBC-V-->> Error: The current user do not have write permissions on the output path: \n'));
                    disp(app.gener_OutputfolderEditField.Value);
                    disp('Please check the folder permission.');
                    checked = false;
                end
            end
            % if(isequal(app.sens_FieldtriptemplateDropDown.Value,'-Select a layout-'))
            %     msgbox('The fieldtrip layout in sensor level tab must be a correct template.','Info');
            %     checked = false;
            % end
            if(isequal(app.delta_runSwitch.Value,'No') && isequal(app.theta_runSwitch.Value,'No')...
                    && isequal(app.alpha_runSwitch.Value,'No') && isequal(app.beta_runSwitch.Value,'No')...
                    && isequal(app.gamma1_runSwitch.Value,'No') && isequal(app.gamma2_runSwitch.Value,'No'))
            msgbox('Please, you have to select at least a frequency to run.','Info');
            checked = false;
            end
        end % check_user_params
        
        % Setting properties files with user interface configuration
        function update_property_files(app)                       
            % Updating general params            
            if(isequal(app.gener_RuninGPUSwitch.Value,'Yes')); app.bcv_properties.general_params.use_gpu.value = true;
            else; app.bcv_properties.general_params.use_gpu.value = false; end
            
            app.bcv_properties.general_params.params.bcv_workspace.Dataset_name            = app.gener_DatasetNameEditField.Value;
            app.bcv_properties.general_params.params.bcv_workspace.Dataset_description     = app.gener_DescriptionTextArea.Value;
            app.bcv_properties.general_params.params.bcv_workspace.BCV_input_dir           = app.gener_InputfolderEditField.Value;
            app.bcv_properties.general_params.params.bcv_workspace.BCV_work_dir            = app.gener_OutputfolderEditField.Value;
            app.bcv_properties.general_params.params.analysis_level.value                  = app.gener_AnalysislevelDropDown.Value;
            app.bcv_properties.general_params.params.analysis_level.reset_all              = true;
            if(isequal(app.gener_RunbytrialsSwitch.Value,'Yes')); app.bcv_properties.general_params.params.run_by_trial.value = true;
            else; app.bcv_properties.general_params.params.run_by_trial.value = false;end
            if(isequal(app.gener_RunbyfrequencybinSwitch.Value,'Yes')); app.bcv_properties.general_params.params.run_frequency_bin.value = true;
            else; app.bcv_properties.general_params.params.run_frequency_bin.value = false; end
            if(isequal(app.gener_GetsystemresponseSwitch.Value,'Yes')); app.bcv_properties.general_params.params.gener_GetsystemresponseSwitch.value = true;
            else; app.bcv_properties.general_params.params.gener_GetsystemresponseSwitch.value = false; end
            
            % Updating sensor level params
            app.bcv_properties.sensor_params.params.method.value                              = app.sens_FFTMethodDropDown.Value;
            % app.bcv_properties.sensor_params.params.fieldtrip.layout.value                  = app.sens_FieldtriptemplateDropDown.Value;

            app.bcv_properties.sensor_params.params.frequencies(1).f_start                    = app.delta_startEditField.Value;
            app.bcv_properties.sensor_params.params.frequencies(1).f_end                      = app.delta_endEditField.Value;
            if(isequal(app.delta_runSwitch.Value,'Yes')); app.bcv_properties.sensor_params.params.frequencies(1).run = true;
            else; app.bcv_properties.sensor_params.params.frequencies(1).run = false; end
            app.bcv_properties.sensor_params.params.frequencies(2).f_start                    = app.theta_startEditField.Value;
            app.bcv_properties.sensor_params.params.frequencies(2).f_end                      = app.theta_endEditField.Value;
            if(isequal(app.theta_runSwitch.Value,'Yes')); app.bcv_properties.sensor_params.params.frequencies(2).run = true;
            else; app.bcv_properties.sensor_params.params.frequencies(2).run = false; end
            app.bcv_properties.sensor_params.params.frequencies(3).f_start                    = app.alpha_startEditField.Value;
            app.bcv_properties.sensor_params.params.frequencies(3).f_end                      = app.alpha_endEditField.Value;
            if(isequal(app.alpha_runSwitch.Value,'Yes')); app.bcv_properties.sensor_params.params.frequencies(3).run = true;
            else; app.bcv_properties.sensor_params.params.frequencies(3).run = false; end
            app.bcv_properties.sensor_params.params.frequencies(4).f_start                    = app.beta_startEditField.Value;
            app.bcv_properties.sensor_params.params.frequencies(4).f_end                      = app.beta_endEditField.Value;
            if(isequal(app.beta_runSwitch.Value,'Yes')); app.bcv_properties.sensor_params.params.frequencies(4).run = true;
            else; app.bcv_properties.sensor_params.params.frequencies(4).run = false; end
            app.bcv_properties.sensor_params.params.frequencies(5).f_start                    = app.gamma1_startEditField.Value;
            app.bcv_properties.sensor_params.params.frequencies(5).f_end                      = app.gamma1_endEditField.Value;
            if(isequal(app.gamma1_runSwitch.Value,'Yes')); app.bcv_properties.sensor_params.params.frequencies(5).run = true;
            else; app.bcv_properties.sensor_params.params.frequencies(5).run = false; end
            app.bcv_properties.sensor_params.params.frequencies(6).f_start                    = app.gamma2_startEditField.Value;
            app.bcv_properties.sensor_params.params.frequencies(6).f_end                      = app.gamma2_endEditField.Value;
            if(isequal(app.gamma2_runSwitch.Value,'Yes')); app.bcv_properties.sensor_params.params.frequencies(6).run = true;
            else; app.bcv_properties.sensor_params.params.frequencies(6).run = false; end
            app.bcv_properties.sensor_params.params.freq_resol.value                          = app.spect_freq_resolEditField.Value;
            app.bcv_properties.sensor_params.params.samp_freq.value                           = app.spect_sample_freqEditField.Value;
            app.bcv_properties.sensor_params.params.max_freq.value                            = app.spect_FreqmaxEditField.Value;
            app.bcv_properties.sensor_params.params.freq_gfiltvar.value                       = app.spect_freq_gfiltvarEditField.Value;
            app.bcv_properties.sensor_params.params.win_order.value                           = app.spect_win_orderSlider.Value;
            
            % Updating activation level params
            app.bcv_properties.activation_params.params.methods{1, 1}.sssblpp.run               = app.CheckBox_sssblpp.Value;
            app.bcv_properties.activation_params.params.methods{1, 1}.sssblpp.sssblpp_th.value  = app.sssblpp_thEditField.Value;
            app.bcv_properties.activation_params.params.methods{2, 1}.eloreta.run               = app.CheckBox_eloreta.Value;
            app.bcv_properties.activation_params.params.methods{2, 1}.eloreta.gamma1.value      = app.eloreta_gamma1EditField.Value;
            app.bcv_properties.activation_params.params.methods{2, 1}.eloreta.gamma2.value      = app.eloreta_gamma2EditField.Value;
            app.bcv_properties.activation_params.params.methods{2, 1}.eloreta.delta_gamma.value = app.eloreta_deltagammaEditField.Value;
            app.bcv_properties.activation_params.params.methods{2, 1}.eloreta.eloreta_th.value  = app.eloreta_thEditField.Value;
            app.bcv_properties.activation_params.params.methods{3, 1}.lcmv.run                  = app.CheckBox_lcmv.Value;
            app.bcv_properties.activation_params.params.methods{3, 1}.lcmv.gamma1.value         = app.lcmv_gamma1EditField.Value;
            app.bcv_properties.activation_params.params.methods{3, 1}.lcmv.gamma2.value         = app.lcmv_gamma2EditField.Value;
            app.bcv_properties.activation_params.params.methods{3, 1}.lcmv.delta_gamma.value    = app.lcmv_deltagammaEditField.Value;
            app.bcv_properties.activation_params.params.methods{3, 1}.lcmv.lcmv_th.value        = app.lcmv_thEditField.Value;
            app.bcv_properties.activation_params.params.IsCurv.value                            = app.activ_IsCurvCheckBox.Value;
            app.bcv_properties.activation_params.params.IsParcel.value                          = app.activ_IsParcelCheckBox.Value;
            app.bcv_properties.activation_params.params.IsNeigh.value                           = app.activ_IsNeighCheckBox.Value;
            app.bcv_properties.activation_params.params.IsField.value                           = app.activ_fieldSlider.Value;
            app.bcv_properties.activation_params.params.aGiri.value                             = app.activ_aGiriEditField.Value;
            app.bcv_properties.activation_params.params.aSulc.value                             = app.activ_aSulcEditField.Value;
            app.bcv_properties.activation_params.params.bGiri.value                             = app.activ_bGiriEditField.Value;
            app.bcv_properties.activation_params.params.bSulc.value                             = app.activ_bSulcEditField.Value;
            app.bcv_properties.activation_params.params.regLaplacian.value                      = app.activ_regLaplacianEditField.Value;
            
            % Updating connectivity level params
            app.bcv_properties.connectivity_params.params.methods{1, 1}.higgs.run                  = app.CheckBox_higgs.Value;
            app.bcv_properties.connectivity_params.params.methods{1, 1}.higgs.higgs_th.value       = app.higgs_thEditField.Value;
            app.bcv_properties.connectivity_params.params.methods{2, 1}.hg_lasso.run               = app.CheckBox_hg_lasso.Value;
            app.bcv_properties.connectivity_params.params.methods{2, 1}.hg_lasso.hg_lasso_th.value = app.hg_lasso_thEditField.Value;
            app.bcv_properties.connectivity_params.params.IsCurv.value                             = app.conn_IsCurvCheckBox.Value;
            app.bcv_properties.connectivity_params.params.IsNeigh.value                            = app.conn_IsNeighCheckBox.Value;
            app.bcv_properties.connectivity_params.params.IsField.value                            = app.conn_fieldSlider.Value;
            app.bcv_properties.connectivity_params.params.aGiri.value                              = app.conn_aGiriEditField.Value;
            app.bcv_properties.connectivity_params.params.aSulc.value                              = app.conn_aSulcEditField.Value;
            app.bcv_properties.connectivity_params.params.bGiri.value                              = app.conn_bGiriEditField.Value;
            app.bcv_properties.connectivity_params.params.bSulc.value                              = app.conn_bSulcEditField.Value;
            app.bcv_properties.connectivity_params.params.axi.value                                = app.conn_axiEditField.Value;
            app.bcv_properties.connectivity_params.params.maxiter_outer.value                      = app.conn_MaxiterouterEditField.Value;
            app.bcv_properties.connectivity_params.params.maxiter_inner.value                      = app.conn_MaxiterinnerEditField.Value;
            app.bcv_properties.connectivity_params.params.ntry.value                               = app.conn_ntrySlider.Value;
            app.bcv_properties.connectivity_params.params.prew.value                               = app.conn_prewarmingCheckBox.Value;
            app.bcv_properties.connectivity_params.params.penalty.value                            = app.conn_penaltySlider.Value;
            app.bcv_properties.connectivity_params.params.rth1.value                               = app.conn_rth1EditField.Value;
            app.bcv_properties.connectivity_params.params.rth2.value                               = app.conn_rth2EditField.Value;
            app.bcv_properties.connectivity_params.params.eigreg.value                             = app.conn_eigregEditField.Value;
            app.bcv_properties.connectivity_params.params.regLaplacian.value                       = app.conn_regLaplacianEditField.Value;
                       
            % saving property files
            saveJSON(app.bcv_properties.general_params,'bcv_properties/general_params.json');
            h = matlab.desktop.editor.openDocument(fullfile(pwd,'bcv_properties/general_params.json'));
            h.smartIndentContents
            h.save
            h.close
            saveJSON(app.bcv_properties.sensor_params,'bcv_properties/sensor_params.json');
            h = matlab.desktop.editor.openDocument(fullfile(pwd,'bcv_properties/sensor_params.json'));
            h.smartIndentContents
            h.save
            h.close
            saveJSON(app.bcv_properties.activation_params,'bcv_properties/activation_params.json');
            h = matlab.desktop.editor.openDocument(fullfile(pwd,'bcv_properties/activation_params.json'));
            h.smartIndentContents
            h.save
            h.close
            saveJSON(app.bcv_properties.connectivity_params,'bcv_properties/connectivity_params.json');
            h = matlab.desktop.editor.openDocument(fullfile(pwd,'bcv_properties/connectivity_params.json'));
            h.smartIndentContents
            h.save
            h.close
            
        end % update_property_files
        
        function load_properties(app)
            app.bcv_properties          = get_properties();            
            app.app_params.generals     = app.bcv_properties.generals;
        end
        
        function set_component_values(app)
            % Loadding App properties
            
            toolbox_name = app.app_params.generals.name;
            version_date = app.app_params.generals.version_date;
            version      = app.app_params.generals.version;
            app.BCVARETAToolboxv10UIFigure.Name = strcat(toolbox_name," ", version, " (", version_date,")");
            
            % Loadding General properties   
            Dataset_name = app.bcv_properties.general_params.params.bcv_workspace.Dataset_name;
            app.gener_DatasetNameEditField.Value = Dataset_name;
            Dataset_description = app.bcv_properties.general_params.params.bcv_workspace.Dataset_name;            
            app.gener_DescriptionTextArea.Value = Dataset_description;
            BCV_input_dir = app.bcv_properties.general_params.params.bcv_workspace.BCV_input_dir;
            app.gener_InputfolderEditField.Value = BCV_input_dir;
            BCV_work_dir = app.bcv_properties.general_params.params.bcv_workspace.BCV_work_dir;
            app.gener_OutputfolderEditField.Value = BCV_work_dir;
            use_gpu = app.bcv_properties.general_params.params.use_gpu.value;
            if use_gpu; app.gener_RuninGPUSwitch.Value = 'Yes'; else; app.gener_RuninGPUSwitch.Value = 'No'; end
            try
                count = gpuDeviceCount("available");
                if(isequal(count,0)); app.gener_RuninGPUSwitch.Value = 'No'; end
                if (isequal(app.gener_RuninGPUSwitch.Value,'Yes')); app.gener_RuninGPUSwitch.FontColor = 'b'; app.gener_RuninGPUSwitch.FontWeight = 'bold'; end
            catch
                app.gener_RuninGPUSwitch.Value = 'No';
            end
            
            run_by_trials = app.bcv_properties.general_params.params.run_by_trial.value;
            if run_by_trials; app.gener_RunbytrialsSwitch.Value = 'Yes'; else; app.gener_RunbytrialsSwitch.Value = 'No'; end
            if run_by_trials; app.gener_RunbytrialsSwitch.FontColor = 'b'; app.gener_RunbytrialsSwitch.FontWeight = 'bold'; end
            run_by_bin = app.bcv_properties.general_params.params.run_frequency_bin.value;
            if run_by_bin; app.gener_RunbyfrequencybinSwitch.Value = 'Yes'; else; app.gener_RunbyfrequencybinSwitch.Value = 'No'; end
            if run_by_bin; app.gener_RunbyfrequencybinSwitch.FontColor = 'b'; app.gener_RunbyfrequencybinSwitch.FontWeight = 'bold'; end
            system_response = app.bcv_properties.general_params.params.system_response.value;
            if system_response; app.gener_GetsystemresponseSwitch.Value = 'Yes'; else; app.gener_GetsystemresponseSwitch.Value = 'No'; end
            if system_response; app.gener_GetsystemresponseSwitch.FontColor = 'b'; app.gener_GetsystemresponseSwitch.FontWeight = 'bold'; end
            
            %% Loadding Sensor properties
            freq_resol = app.bcv_properties.sensor_params.params.freq_resol.value;
            app.spect_freq_resolEditField.Value = freq_resol;
            samp_freq = app.bcv_properties.sensor_params.params.samp_freq.value;
            app.spect_sample_freqEditField.Value = samp_freq;
            freq_gfiltvar = app.bcv_properties.sensor_params.params.freq_gfiltvar.value;
            app.spect_freq_gfiltvarEditField.Value = freq_gfiltvar;
            max_freq = app.bcv_properties.sensor_params.params.max_freq.value; 
            app.spect_FreqmaxEditField.Value = max_freq;
            win_order = app.bcv_properties.sensor_params.params.win_order.value;
            app.spect_win_orderSlider.Value = win_order;

            frequencies     = app.bcv_properties.sensor_params.params.frequencies;
            for i=1:length(frequencies)
                freq = frequencies(i);
                switch lower(freq.name)
                    case 'delta'
                        if freq.run;  app.delta_runSwitch.Value = 'Yes'; else; app.delta_runSwitch.Value = 'No'; end
                        if freq.run;  app.delta_runSwitch.FontColor = 'b'; app.delta_runSwitch.FontWeight = 'bold'; end
                        if freq.run;  app.DeltaLabel.FontColor = 'b'; app.DeltaLabel.FontWeight = 'bold'; end
                        app.delta_startEditField.Value = freq.f_start;
                        app.delta_endEditField.Value = freq.f_end;
                    case 'theta'
                        if freq.run;  app.theta_runSwitch.Value = 'Yes'; else; app.theta_runSwitch.Value = 'No'; end
                        if freq.run;  app.theta_runSwitch.FontColor = 'b'; app.theta_runSwitch.FontWeight = 'bold'; end
                        if freq.run;  app.ThetaLabel.FontColor = 'b'; app.ThetaLabel.FontWeight = 'bold'; end
                        app.theta_startEditField.Value = freq.f_start;
                        app.theta_endEditField.Value = freq.f_end;
                    case 'alpha'
                        if freq.run;  app.alpha_runSwitch.Value = 'Yes'; else; app.alpha_runSwitch.Value = 'No'; end
                        if freq.run;  app.alpha_runSwitch.FontColor = 'b'; app.alpha_runSwitch.FontWeight = 'bold'; end
                        if freq.run;  app.AlphaLabel.FontColor = 'b'; app.AlphaLabel.FontWeight = 'bold'; end
                        app.alpha_startEditField.Value = freq.f_start;
                        app.alpha_endEditField.Value = freq.f_end;
                    case 'beta'
                        if freq.run;  app.beta_runSwitch.Value = 'Yes'; else; app.beta_runSwitch.Value = 'No'; end
                        if freq.run;  app.beta_runSwitch.FontColor = 'b'; app.beta_runSwitch.FontWeight = 'bold'; end
                        if freq.run;  app.BetaLabel.FontColor = 'b'; app.BetaLabel.FontWeight = 'bold'; end
                        app.beta_startEditField.Value = freq.f_start;
                        app.beta_endEditField.Value = freq.f_end;
                    case 'gamma1'
                        if freq.run;  app.gamma1_runSwitch.Value = 'Yes'; else; app.gamma1_runSwitch.Value = 'No'; end
                        if freq.run;  app.gamma1_runSwitch.FontColor = 'b'; app.gamma1_runSwitch.FontWeight = 'bold'; end
                        if freq.run;  app.Gamma1Label.FontColor = 'b'; app.Gamma1Label.FontWeight = 'bold'; end
                        app.gamma1_startEditField.Value = freq.f_start;
                        app.gamma1_endEditField.Value = freq.f_end;
                    case 'gamma2'
                        if freq.run;  app.gamma2_runSwitch.Value = 'Yes'; else; app.gamma2_runSwitch.Value = 'No'; end
                        if freq.run;  app.gamma2_runSwitch.FontColor = 'b'; app.gamma2_runSwitch.FontWeight = 'bold'; end
                        if freq.run;  app.Gamma2Label.FontColor = 'b'; app.Gamma2Label.FontWeight = 'bold'; end
                        app.gamma2_startEditField.Value = freq.f_start;
                        app.gamma2_endEditField.Value = freq.f_end;
                end
            end
            

            % ft_template_dir = fullfile("external/fieldtrip/template/layout");
            % layouts = dir(ft_template_dir);
            % layouts(ismember( {layouts.name}, {'.', '..'})) = [];  %remove . and ..
            % app.sens_FieldtriptemplateDropDown.Items = {layouts.name};
            % app.sens_FieldtriptemplateDropDown.Items{end+1} = '-Select a layout-';
            % app.sens_FieldtriptemplateDropDown.Value = '-Select a layout-';
            % 
            % fieldtrip_layout = app.sensor_params.params.fieldtrip.layout.value;
            % if(~isempty(find(contains( {layouts.name}, fieldtrip_layout),1)))
            %     app.sens_FieldtriptemplateDropDown.Value = fieldtrip_layout;
            % elseif(~isequal(lower(fieldtrip_layout),'none') || isequal(lower(fieldtrip_layout),'default') || ~isempty(fieldtrip_layout))
            %     app.sens_FieldtriptemplateDropDown.Value = '-Select a layout-';
            % else
            %     msgbox({'The selected layout is not contained in the fildtrip layouts.',...
            %         fieldtrip_layout,...
            %         'Please configure it in the bcv_properties/sensor_level.json file.'},'Info');
            % end
            
            %% Loadding Activation properties
            
            activ_methods = app.bcv_properties.activation_params.params.methods;
            for i=1:length(activ_methods)
                method = activ_methods{i};
                method_name = fieldnames(method);
                method_name = method_name{1};
                switch method_name
                    case 'sssblpp'
                        run = method.(method_name).run;
                        app.CheckBox_sssblpp.Value = run;
                        threshold = method.(method_name).sssblpp_th.value;
                        app.sssblpp_thEditField.Value = threshold;
                    case 'eloreta'
                        run = method.(method_name).run;
                        app.CheckBox_eloreta.Value = run;
                        threshold = method.(method_name).eloreta_th.value;
                        app.eloreta_thEditField.Value = threshold;
                        gamma1 = method.(method_name).gamma1.value;
                        app.eloreta_gamma1EditField.Value = gamma1;
                        gamma2 = method.(method_name).gamma2.value;
                        app.eloreta_gamma2EditField.Value = gamma2;
                        delta_gamma = method.(method_name).delta_gamma.value;
                        app.eloreta_deltagammaEditField.Value = delta_gamma;
                    case 'lcmv'
                        run = method.(method_name).run;
                        app.CheckBox_lcmv.Value = run;
                        threshold = method.(method_name).lcmv_th.value;
                        app.lcmv_thEditField.Value = threshold;
                        gamma1 = method.(method_name).gamma1.value;
                        app.lcmv_gamma1EditField.Value = gamma1;
                        gamma2 = method.(method_name).gamma2.value;
                        app.lcmv_gamma2EditField.Value = gamma2;
                        delta_gamma = method.(method_name).delta_gamma.value;
                        app.lcmv_deltagammaEditField.Value = delta_gamma;
                end
            end
            activ_IsCurv = app.bcv_properties.activation_params.params.IsCurv.value;
            app.activ_IsCurvCheckBox.Value = activ_IsCurv;
            activ_IsParcel = app.bcv_properties.activation_params.params.IsParcel.value;
            app.activ_IsParcelCheckBox.Value = activ_IsParcel;
            activ_IsNeigh = app.bcv_properties.activation_params.params.IsNeigh.value;
            app.activ_IsNeighCheckBox.Value = activ_IsNeigh;
            activ_aGiri = app.bcv_properties.activation_params.params.aGiri.value;
            app.activ_aGiriEditField.Value = activ_aGiri;
            activ_bGiri = app.bcv_properties.activation_params.params.bGiri.value;
            app.activ_bGiriEditField.Value = activ_bGiri;
            activ_aSulc = app.bcv_properties.activation_params.params.aSulc.value;
            app.activ_aSulcEditField.Value = activ_aSulc;
            activ_bSulc = app.bcv_properties.activation_params.params.bSulc.value;
            app.activ_bSulcEditField.Value = activ_bSulc;
            activ_Field = app.bcv_properties.activation_params.params.IsField.value;
            app.activ_fieldSlider.Value = activ_Field;
            activ_regLaplacian = app.bcv_properties.activation_params.params.regLaplacian.value;
            app.activ_regLaplacianEditField.Value = activ_regLaplacian;
            
            %% Loadding Connectivity properties
            conn_methods = app.bcv_properties.connectivity_params.params.methods;
            for i=1:length(conn_methods)
                method = conn_methods{i};
                method_name = fieldnames(method);
                method_name = method_name{1};
                switch method_name
                    case 'higgs'
                        run = method.(method_name).run;
                        app.CheckBox_higgs.Value = run;
                        threshold = method.(method_name).higgs_th.value;
                        app.higgs_thEditField.Value = threshold;
                    case 'hg_lasso'
                        run = method.(method_name).run;
                        app.CheckBox_hg_lasso.Value = run;
                        threshold = method.(method_name).hg_lasso_th.value;
                        app.hg_lasso_thEditField.Value = threshold;
                end
            end
            conn_IsCurv = app.bcv_properties.connectivity_params.params.IsCurv.value;
            app.conn_IsCurvCheckBox.Value = conn_IsCurv;
            conn_IsNeigh = app.bcv_properties.connectivity_params.params.IsNeigh.value;
            app.conn_IsNeighCheckBox.Value = conn_IsNeigh;
            conn_prew = app.bcv_properties.connectivity_params.params.prew.value;
            app.conn_prewarmingCheckBox.Value = conn_prew;
            conn_aGiri = app.bcv_properties.connectivity_params.params.aGiri.value;
            app.conn_aGiriEditField.Value = conn_aGiri;
            conn_bGiri = app.bcv_properties.connectivity_params.params.bGiri.value;
            app.conn_bGiriEditField.Value = conn_bGiri;
            conn_aSulc = app.bcv_properties.connectivity_params.params.aSulc.value;
            app.conn_aSulcEditField.Value = conn_aSulc;
            conn_bSulc = app.bcv_properties.connectivity_params.params.bSulc.value;
            app.conn_bSulcEditField.Value = conn_bSulc;
            conn_outer = app.bcv_properties.connectivity_params.params.maxiter_outer.value;
            app.conn_MaxiterouterEditField.Value = conn_outer;
            conn_inner = app.bcv_properties.connectivity_params.params.maxiter_inner.value;
            app.conn_MaxiterinnerEditField.Value = conn_inner;
            conn_axi = app.bcv_properties.connectivity_params.params.axi.value;
            app.conn_axiEditField.Value = conn_axi;
            conn_field = app.bcv_properties.connectivity_params.params.IsField.value;
            app.conn_fieldSlider.Value = conn_field;
            conn_penalty = app.bcv_properties.connectivity_params.params.penalty.value;
            app.conn_penaltySlider.Value = conn_penalty;
            conn_rth1 = app.bcv_properties.connectivity_params.params.rth1.value;
            app.conn_rth1EditField.Value = conn_rth1;
            conn_rth2 = app.bcv_properties.connectivity_params.params.rth2.value;
            app.conn_rth2EditField.Value = conn_rth2;
            conn_ntry = app.bcv_properties.connectivity_params.params.ntry.value;
            app.conn_ntrySlider.Value = conn_ntry;
            conn_eigreg = app.bcv_properties.connectivity_params.params.eigreg.value;
            app.conn_eigregEditField.Value = conn_eigreg;
            conn_regLaplacian = app.bcv_properties.connectivity_params.params.regLaplacian.value;
            app.conn_regLaplacianEditField.Value = conn_regLaplacian; 
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            %% Creating BC-VARETA files
            homedir = char(java.lang.System.getProperty('user.home'));
            app.BCVdir  = fullfile(homedir,"BC_VARETA");
            if(~isfolder(app.BCVdir))
                mkdir(app.BCVdir);
            end
            app.Datasets_file = fullfile(app.BCVdir,"Datasets.json");
            if(isfile(app.Datasets_file))
                TempDatasets = jsondecode(fileread(app.Datasets_file));
                if(~isempty(TempDatasets))
                    opts.Interpreter = 'tex';
                    % Include the desired Default answer
                    opts.Default = 'Yes';
                    % Use the TeX interpreter to format the question
                    quest = {strcat("BC-VARETA db have saved Datasets."), ...
                        'Would you like import this the datasets?'};
                    answer = questdlg(quest,'Apply changes',...
                        'Load','Cancel',opts);
                    % Handle response
                    switch answer
                        case 'Load'
                            app.Datasets = TempDatasets;
                            [icondata,iconcmap] = imread("guide/images/success.png");
                            msg = msgbox({'The Datasets loaded successfully!!!'},'Success',"custom",icondata,iconcmap);
                            movegui(msg,'north');
                            disp("-->> Datasets loaded");                        
                        case 'Cancel'
                            app.Datasets = struct;
                            saveJSON(app.Datasets,app.Datasets_file);
                    end
                end
            else
                app.Datasets = struct([]);
                saveJSON(app.Datasets,app.Datasets_file);
            end
            %% Load properties from files
            load_properties(app);  
            %% Set components from files
            set_component_values(app); 
            %             try
            %                 jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
            %                 jCmdWin = jDesktop.getClient('Command Window');
            %                 jTextArea = jCmdWin.getComponent(0).getViewport.getView;
            %                 set(jTextArea,'CaretUpdateCallback',@app.setPromptFcn)
            %             catch
            %                 warndlg('fatal error');
            %             end
        end

        % Menu selected function: ExitMenu
        function ExitMenuSelected(app, event)
            delete(app);
        end

        % Button pushed function: gener_SelectInputButton
        function gener_SelectInputButtonPushed(app, event)
            try
                folder = uigetdir("Select the input data folder");
                subjects = dir(fullfile(folder,'**','subject.mat'));
                if(~isempty(subjects))
                    app.gener_InputfolderEditField.Value = folder;
                else
                    msg = msgbox({'The selected folder do not contains any BC-VARETA structure.',...
                        ' Please select a correct folder.'},'Info');
                    movegui(msg,'north');
                end
            catch
            end
        end

        % Button pushed function: gener_SelectOutputButton
        function gener_SelectOutputButtonPushed(app, event)
            try
                folder = uigetdir("Select the Output folder");
                [~,values] = fileattrib(folder);
                if(values.UserWrite)
                    app.gener_OutputfolderEditField.Value = folder;
                else
                    msg = msgbox({'The current user do not have write permissions on the selected forder.',...
                        ' Please check the folder permission or select another output folder.'},'Info');
                    movegui(msg,'north');
                end
            catch
            end
        end

        % Value changed function: gener_AnalysislevelDropDown
        function gener_AnalysislevelDropDownValueChanged(app, event)
            analysis_level = app.gener_AnalysislevelDropDown.Value;
            switch analysis_level
                case 'all'
                    app.SensorTab.Parent = [];
                    app.ActivationTab.Parent = [];
                    app.ConnectivityTab.Parent = [];                    
                    app.SensorTab.Parent = app.TabGroup;
                    app.ActivationTab.Parent = app.TabGroup;
                    app.ConnectivityTab.Parent = app.TabGroup;                   
                case '1'
                    app.SensorTab.Parent = [];
                    app.ActivationTab.Parent = [];
                    app.ConnectivityTab.Parent = [];                    
                    app.SensorTab.Parent = app.TabGroup;                   
                case '2'
                    app.SensorTab.Parent = [];
                    app.ActivationTab.Parent = [];
                    app.ConnectivityTab.Parent = [];                   
                    app.ActivationTab.Parent = app.TabGroup;
                case '3'
                    app.SensorTab.Parent = [];
                    app.ActivationTab.Parent = [];
                    app.ConnectivityTab.Parent = [];                    
                    app.ConnectivityTab.Parent = app.TabGroup;
                case '12'
                    app.SensorTab.Parent = [];
                    app.ActivationTab.Parent = [];
                    app.ConnectivityTab.Parent = [];
                    app.SensorTab.Parent = app.TabGroup;
                    app.ActivationTab.Parent = app.TabGroup;
                case '23'
                    app.SensorTab.Parent = [];
                    app.ActivationTab.Parent = [];
                    app.ConnectivityTab.Parent = [];
                    app.ActivationTab.Parent = app.TabGroup;
                    app.ConnectivityTab.Parent = app.TabGroup;
            end
        end

        % Button pushed function: RunButton
        function RunButtonPushed(app, event)
            checked = check_user_params(app);
            if checked
                update_property_files(app); 
                opts.Interpreter = 'tex';
                % Include the desired Default answer
                opts.Default = 'Yes';
                % Use the TeX interpreter to format the question
                quest = {strcat("You will run: ",app.BCVARETAToolboxv10UIFigure.Name,"."), ...
                    'Would you like continue the process?'};
                answer = questdlg(quest,'Start process',...
                    'Yes','No',opts);
                % Handle response
                switch answer
                    case 'Yes'
                        % f = dialog('Position',[300 300 250 80]);
                        % movegui(f,'center')
                        % iconsClassName = 'com.mathworks.widgets.BusyAffordance$AffordanceSize';
                        % iconsSizeEnums = javaMethod('values',iconsClassName);
                        % SIZE_32x32 = iconsSizeEnums(2);  % (1) = 16x16,  (2) = 32x32
                        % jObj = com.mathworks.widgets.BusyAffordance(SIZE_32x32, 'Processing');  % icon, label
                        % jObj.setPaintsWhenStopped(true);  % default = false
                        % jObj.useWhiteDots(false);         % default = false (true is good for dark backgrounds)
                        % javacomponent(jObj.getComponent, [50,10,150,80], f);
                        % jObj.start;
                        % pause(1);                           
                        BC_VARETA_bash(1,1);
                        % jObj.stop;
                        % jObj.setBusyText('All done!');
                        % disp('All done....');
                        % pause(2);
                        % delete(f);
                        % msgbox({'BC-VARETA Toolbox process completed.'},'Info');
                    case 'No'
                        
                end
                
                
            end
        end

        % Value changed function: gener_RunbytrialsSwitch
        function gener_RunbytrialsSwitchValueChanged(app, event)
            value = app.gener_RunbytrialsSwitch.Value;
            if(isequal(value,'Yes'))
                try
                    subject_files = dir(fullfile(app.gener_InputfolderEditField.Value,'**','subject.mat'));
                    file_test = subject_files(1);
                    load(fullfile(file_test.folder,file_test.name));
                    load(fullfile(file_test.folder,subject_info.meeg_dir));
                    if(~iscell(MEEG.data))
                        msgbox({'The selected input data can not be run by trials.',...
                            'Please check the input dataset or Run BC-VARETA without trials.'},'Info');
                        app.gener_RunbytrialsSwitch.Value = 'No';
                        app.gener_RunbytrialsSwitch.FontColor = 'k';
                        app.gener_RunbytrialsSwitch.FontWeight = 'normal';
                    else
                        app.gener_RunbytrialsSwitch.FontColor = 'b';
                        app.gener_RunbytrialsSwitch.FontWeight = 'bold';
                    end
                catch
                    msg = msgbox('Please select first the input data folder.','Info');
                    movegui(msg,'north');
                    app.gener_RunbytrialsSwitch.Value = 'No';
                end
            else
                app.gener_RunbytrialsSwitch.FontColor = 'k';
                app.gener_RunbytrialsSwitch.FontWeight = 'normal';
            end
        end

        % Value changed function: gener_RuninGPUSwitch
        function gener_RuninGPUSwitchValueChanged(app, event)
            value = app.gener_RuninGPUSwitch.Value;
            if(isequal(value,'Yes'))
                try
                    count = gpuDeviceCount();
                    if(isequal(count,0))
                        msgbox({'BC-VARETA Toolbox con not be run in GPU on this PC.',...
                            'This PC do not have any GPU available.'},'Info');
                        app.gener_RuninGPUSwitch.Value = 'No';
                    else
                        gpuinfo = gpuDevice();
                        opts.Interpreter = 'tex';
                        % Include the desired Default answer
                        opts.Default = 'Yes';
                        % Use the TeX interpreter to format the question
                        quest = {strcat("GPU divice detected: ",gpuinfo.Name,"."), ...
                            strcat("Total memory: ",num2str(gpuinfo.TotalMemory/(1024^3))," GB."), ...
                            strcat("Available memory: ",num2str(gpuinfo.AvailableMemory/(1024^3))," GB."),...
                            strcat("Note: The GPU process was tested with 6GB of memory."),...
                            'Would you like to run BC-VARETA in GPU any way?'};
                        answer = questdlg(quest,'Use GPU',...
                            'Yes','No',opts);
                        % Handle response
                        switch answer
                            case 'Yes'
                                app.gener_RuninGPUSwitch.Value = 'Yes';
                                app.gener_RuninGPUSwitch.FontColor = 'b';
                                app.gener_RuninGPUSwitch.FontWeight = 'bold';
                            case 'No'
                                app.gener_RuninGPUSwitch.Value = 'No';
                                app.gener_RuninGPUSwitch.FontColor = 'k';
                                app.gener_RuninGPUSwitch.FontWeight = 'normal';
                        end
                    end
                catch
                    app.gener_RuninGPUSwitch.FontColor = 'k';
                    app.gener_RuninGPUSwitch.FontWeight = 'normal';
                    msg = msgbox({'BC-VARETA Toolbox con not be run in GPU on this PC.',...
                        'This PC do not have any GPU available.'},'Info');
                    movegui(msg,'north');
                    app.gener_RuninGPUSwitch.Value = 'No';
                end
            else
                app.gener_RuninGPUSwitch.FontColor = 'k';
                app.gener_RuninGPUSwitch.FontWeight = 'normal';
            end
        end

        % Value changed function: delta_endEditField
        function delta_endEditFieldValueChanged(app, event)
            end_value = app.delta_endEditField.Value;
            if(end_value<= app.delta_startEditField.Value)
                msg = msgbox({'The end frequency value must be greater than the start value.',...
                    'Please check the frequency input values.'},'Info');
                movegui(msg,'north');
                app.delta_endEditField.Value = event.PreviousValue;
            end
        end

        % Value changed function: delta_startEditField
        function delta_startEditFieldValueChanged(app, event)
            start_value = app.delta_startEditField.Value;
            if(start_value>= app.delta_endEditField.Value)
                msg = msgbox({'The start frequency value must be less than the end value.',...
                    'Please check the frequency input values.'},'Info');
                movegui(msg,'north');
                app.delta_startEditField.Value = event.PreviousValue;
            end
        end

        % Value changed function: theta_startEditField
        function theta_startEditFieldValueChanged(app, event)
            start_value = app.theta_startEditField.Value;
            if(start_value>= app.theta_endEditField.Value)
                msg = msgbox({'The start frequency value must be less than the end value.',...
                    'Please check the frequency input values.'},'Info');
                movegui(msg,'north');
                app.theta_startEditField.Value = event.PreviousValue;
            end
        end

        % Value changed function: theta_endEditField
        function theta_endEditFieldValueChanged(app, event)
            end_value = app.theta_endEditField.Value;
            if(end_value<= app.theta_startEditField.Value)
                msg = msgbox({'The end frequency value must be greater than the start value.',...
                    'Please check the frequency input values.'},'Info');
                movegui(msg,'north');
                app.theta_endEditField.Value = event.PreviousValue;
            end
        end

        % Value changed function: alpha_startEditField
        function alpha_startEditFieldValueChanged(app, event)
            start_value = app.alpha_startEditField.Value;
            if(start_value>= app.alpha_endEditField.Value)
                msg = msgbox({'The start frequency value must be less than the end value.',...
                    'Please check the frequency input values.'},'Info');
                movegui(msg,'north');
                app.alpha_startEditField.Value = event.PreviousValue;
            end
        end

        % Value changed function: alpha_endEditField
        function alpha_endEditFieldValueChanged(app, event)
            end_value = app.alpha_endEditField.Value;
            if(end_value<= app.alpha_startEditField.Value)
                msg = msgbox({'The end frequency value must be greater than the start value.',...
                    'Please check the frequency input values.'},'Info');
                movegui(msg,'north');
                app.alpha_endEditField.Value = event.PreviousValue;
            end
        end

        % Value changed function: beta_startEditField
        function beta_startEditFieldValueChanged(app, event)
            start_value = app.beta_startEditField.Value;
            if(start_value>= app.beta_endEditField.Value)
                msg = msgbox({'The start frequency value must be less than the end value.',...
                    'Please check the frequency input values.'},'Info');
                movegui(msg,'north');
                app.beta_startEditField.Value = event.PreviousValue;
            end
        end

        % Value changed function: beta_endEditField
        function beta_endEditFieldValueChanged(app, event)
            end_value = app.beta_endEditField.Value;
            if(end_value<= app.beta_startEditField.Value)
                msg = msgbox({'The end frequency value must be greater than the start value.',...
                    'Please check the frequency input values.'},'Info');
                movegui(msg,'north');
                app.beta_endEditField.Value = event.PreviousValue;
            end
        end

        % Value changed function: gamma1_startEditField
        function gamma1_startEditFieldValueChanged(app, event)
            start_value = app.gamma1_startEditField.Value;
            if(start_value>= app.gamma1_endEditField.Value)
                msg = msgbox({'The start frequency value must be less than the end value.',...
                    'Please check the frequency input values.'},'Info');
                movegui(msg,'north');
                app.gamma1_startEditField.Value = event.PreviousValue;
            end
        end

        % Value changed function: gamma1_endEditField
        function gamma1_endEditFieldValueChanged(app, event)
            end_value = app.gamma1_endEditField.Value;
            if(end_value<= app.gamma1_startEditField.Value)
                msg = msgbox({'The end frequency value must be greater than the start value.',...
                    'Please check the frequency input values.'},'Info');
                movegui(msg,'north');
                app.gamma1_endEditField.Value = event.PreviousValue;
            end
        end

        % Value changed function: gamma2_startEditField
        function gamma2_startEditFieldValueChanged(app, event)
            start_value = app.gamma2_startEditField.Value;
            if(start_value>= app.gamma2_endEditField.Value)
                msg = msgbox({'The start frequency value must be less than the end value.',...
                    'Please check the frequency input values.'},'Info');
                movegui(msg,'north');
                app.gamma2_startEditField.Value = event.PreviousValue;
            end
        end

        % Value changed function: gamma2_endEditField
        function gamma2_endEditFieldValueChanged(app, event)
            end_value = app.gamma2_endEditField.Value;
            if(end_value<= app.gamma2_startEditField.Value)
                msgbox({'The end frequency value must be greater than the start value.',...
                    'Please check the frequency input values.'},'Info');
                app.gamma2_endEditField.Value = event.PreviousValue;
            end
        end

        % Value changed function: delta_runSwitch
        function delta_runSwitchValueChanged(app, event)
            value = app.delta_runSwitch.Value;
            if(isequal(value,'Yes'))
                app.delta_runSwitch.FontColor = 'b';
                app.delta_runSwitch.FontWeight = 'bold';
                app.DeltaLabel.FontColor = 'b';
                app.DeltaLabel.FontWeight = 'bold';
            else
                app.delta_runSwitch.FontColor = 'k';
                app.delta_runSwitch.FontWeight = 'normal';
                app.DeltaLabel.FontColor = 'k';
                app.DeltaLabel.FontWeight = 'normal';
            end
        end

        % Value changed function: theta_runSwitch
        function theta_runSwitchValueChanged(app, event)
            value = app.theta_runSwitch.Value;
            if(isequal(value,'Yes'))
                app.theta_runSwitch.FontColor = 'b';
                app.theta_runSwitch.FontWeight = 'bold';
                app.ThetaLabel.FontColor = 'b';
                app.ThetaLabel.FontWeight = 'bold';
            else
                app.theta_runSwitch.FontColor = 'k';
                app.theta_runSwitch.FontWeight = 'normal';
                app.ThetaLabel.FontColor = 'k';
                app.ThetaLabel.FontWeight = 'normal';
            end
        end

        % Value changed function: alpha_runSwitch
        function alpha_runSwitchValueChanged(app, event)
            value = app.alpha_runSwitch.Value;
            if(isequal(value,'Yes'))
                app.alpha_runSwitch.FontColor = 'b';
                app.alpha_runSwitch.FontWeight = 'bold';
                app.AlphaLabel.FontColor = 'b';
                app.AlphaLabel.FontWeight = 'bold';
            else
                app.alpha_runSwitch.FontColor = 'k';
                app.alpha_runSwitch.FontWeight = 'normal';
                app.AlphaLabel.FontColor = 'k';
                app.AlphaLabel.FontWeight = 'normal';
            end
        end

        % Value changed function: beta_runSwitch
        function beta_runSwitchValueChanged(app, event)
            value = app.beta_runSwitch.Value;
            if(isequal(value,'Yes'))
                app.beta_runSwitch.FontColor = 'b';
                app.beta_runSwitch.FontWeight = 'bold';
                app.BetaLabel.FontColor = 'b';
                app.BetaLabel.FontWeight = 'bold';
            else
                app.beta_runSwitch.FontColor = 'k';
                app.beta_runSwitch.FontWeight = 'normal';
                app.BetaLabel.FontColor = 'k';
                app.BetaLabel.FontWeight = 'normal';
            end
        end

        % Value changed function: gamma1_runSwitch
        function gamma1_runSwitchValueChanged(app, event)
            value = app.gamma1_runSwitch.Value;
            if(isequal(value,'Yes'))
                app.gamma1_runSwitch.FontColor = 'b';
                app.gamma1_runSwitch.FontWeight = 'bold';
                app.Gamma1Label.FontColor = 'b';
                app.Gamma1Label.FontWeight = 'bold';
            else
                app.gamma1_runSwitch.FontColor = 'k';
                app.gamma1_runSwitch.FontWeight = 'normal';
                app.Gamma1Label.FontColor = 'k';
                app.Gamma1Label.FontWeight = 'normal';
            end
        end

        % Value changed function: gamma2_runSwitch
        function gamma2_runSwitchValueChanged(app, event)
            value = app.gamma2_runSwitch.Value;
            if(isequal(value,'Yes'))
                app.gamma2_runSwitch.FontColor = 'b';
                app.gamma2_runSwitch.FontWeight = 'bold';
                app.Gamma2Label.FontColor = 'b';
                app.Gamma2Label.FontWeight = 'bold';
            else
                app.gamma2_runSwitch.FontColor = 'k';
                app.gamma2_runSwitch.FontWeight = 'normal';
                app.Gamma2Label.FontColor = 'k';
                app.Gamma2Label.FontWeight = 'normal';
            end
        end

        % Value changed function: gener_RunbyfrequencybinSwitch
        function gener_RunbyfrequencybinSwitchValueChanged(app, event)
            value = app.gener_RunbyfrequencybinSwitch.Value;
            if(isequal(value,'Yes'))
                app.gener_RunbyfrequencybinSwitch.FontColor = 'b';
                app.gener_RunbyfrequencybinSwitch.FontWeight = 'bold';
            else
                app.gener_RunbyfrequencybinSwitch.FontColor = 'k';
                app.gener_RunbyfrequencybinSwitch.FontWeight = 'normal';
            end
        end

        % Value changed function: gener_GetsystemresponseSwitch
        function gener_GetsystemresponseSwitchValueChanged(app, event)
            value = app.gener_GetsystemresponseSwitch.Value;
            if(isequal(value,'Yes'))
                app.gener_GetsystemresponseSwitch.FontColor = 'b';
                app.gener_GetsystemresponseSwitch.FontWeight = 'bold';
            else
                app.gener_GetsystemresponseSwitch.FontColor = 'k';
                app.gener_GetsystemresponseSwitch.FontWeight = 'normal';
            end
        end

        % Button pushed function: CancelButton
        function CancelButtonPushed(app, event)
            delete(app);
        end

        % Callback function: Default_ParamsPushTool
        function Default_ParamsPushToolClicked(app, event)
            opts.Interpreter = 'tex';
            % Include the desired Default answer
            opts.Default = 'Yes';
            % Use the TeX interpreter to format the question
            quest = {strcat("You will restore the default params"), ...
                'Would you like continue?'};
            answer = questdlg(quest,'Apply changes',...
                'Yes','No',opts);
            % Handle response
            switch answer
                case 'Yes'
                    set_default_properties(app);
                case 'No'
            end
        end

        % Value changed function: gener_OutputfolderEditField
        function gener_OutputfolderEditFieldValueChanged(app, event)
            output_path = app.gener_OutputfolderEditField.Value;
            [~,values] = fileattrib(output_path);
            if(~values.UserWrite)
                msg = msgbox({'The current user do not have write permission in this forlder.',...
                    'Please check user permission or chenge the Output path.'},'Info');
                movegui(msg,'north');
            end
        end

        % Menu selected function: DatasetLoadMenu
        function DatasetLoadMenuSelected(app, event)
            folder = uigetdir("Select the Dataset folder");
            dataset_file = fullfile(folder,'BC_VARETA.mat');
            if(isfile(dataset_file))
                BC_VARETA_info = load(dataset_file);
                BC_VARETA_info.Path = folder;
                if(~isempty(app.Datasets) && contains(BC_VARETA_info.Dataset_name,{app.Datasets.Dataset_name}))
                    opts.Interpreter = 'tex';
                    % Include the desired Default answer
                    opts.Default = 'Yes';
                    % Use the TeX interpreter to format the question
                    quest = {strcat("The dataset: ", BC_VARETA_info.Dataset_name, " is already loaded"), ...
                        'Would you like import this the dataset anyway?'};
                    answer = questdlg(quest,'Apply changes',...
                        'Add','Replace','Cancel',opts);
                    % Handle response
                    switch answer
                        case 'Add'
                            app.Datasets(end + 1) = BC_VARETA_info;
                            [icondata,iconcmap] = imread("guide/images/success.png");
                            msg = msgbox({'The Dataset was loaded successfully!!!'},'Success',"custom",icondata,iconcmap);
                            movegui(msg,'north');
                            disp("-->> Dataset loaded");
                            saveJSON(app.Datasets,app.Datasets_file);
                        case 'Replace'
                            index = find(contains(BC_VARETA_info.Dataset_name,{app.Datasets.Name}),1);
                            app.Datasets(index) = BC_VARETA_info;
                            saveJSON(app.Datasets,app.Datasets_file);
                            disp("-->> Dataset replaced");
                            [icondata,iconcmap] = imread("guide/images/success.png");
                            msg = msgbox({'The Dataset was replaced successfully!!!'},'Success',"custom",icondata,iconcmap);
                            movegui(msg,'north');
                        case 'Cancel'
                    end
                else
                    app.Datasets = BC_VARETA_info;
                    [icondata,iconcmap] = imread("guide/images/success.png");
                    msgbox({'The Dataset was loaded successfully!!!'},'Success',"custom",icondata,iconcmap);
                    disp("-->> Dataset loaded");
                    saveJSON(app.Datasets,app.Datasets_file);
                end
            else
                msg = msgbox({'The selected folder do not contains the dataset file.',...
                    ' Please select a correct folder.'},'Error',"error","modal");
                movegui(msg,'north');
            end

        end

        % Menu selected function: TestDataDownloadMenu
        function TestDataDownloadMenuSelected(app, event)
            folder = uigetdir('tittle','Select the Source Folder');
            if(folder==0)
                return;
            end
            
            f = dialog('Position',[300 300 250 80]);
            
            iconsClassName = 'com.mathworks.widgets.BusyAffordance$AffordanceSize';
            iconsSizeEnums = javaMethod('values',iconsClassName);
            SIZE_32x32 = iconsSizeEnums(2);  % (1) = 16x16,  (2) = 32x32
            jObj = com.mathworks.widgets.BusyAffordance(SIZE_32x32, 'Downloading test data...');  % icon, label
            
            jObj.setPaintsWhenStopped(true);  % default = false
            jObj.useWhiteDots(false);         % default = false (true is good for dark backgrounds)
            javacomponent(jObj.getComponent, [50,10,150,80], f);
            jObj.start;
            pause(1);
            
            app_properties = jsondecode(fileread(strcat('app_properties.json')));
            url = app_properties.generals.test_data_url;
            filename = strcat(folder,filesep,'BC_VARETA_test_data.zip');
            matlab.net.http.HTTPOptions.VerifyServerName = false;
            options = weboptions('Timeout',Inf,'RequestMethod','get');
            
            try
                disp('Downloding test data....');
                outfilename = websave(filename,url,options);
            catch
                delete(f);
                errordlg('Download error!!!','Error');
                return;
            end
            try
                disp('Unpacking test data....');
                exampleFiles = unzip(filename,folder);
            catch
                delete(f);
                errordlg('Unpackage error!!!','Error');
                return;
            end
            jObj.stop;
            jObj.setBusyText('All done!');
            disp('All done....');
            pause(2);
            delete(f);
            msgbox('Download complete','Info');
        end

        % Menu selected function: TestDataImportMenu
        function TestDataImportMenuSelected(app, event)
            import_data_structure;
        end

        % Menu selected function: DatasetListMenu
        function DatasetListMenuSelected(app, event)
            if(isempty(app.Datasets))
                msg = msgbox({'There is no Dataset to show.',...
                    'Please import or create a Dataset.'},'Info',"help","modal");
                movegui(msg,'north');
            else
                datasetV                    = Dataset;    
                datasetV.BCVdir             = app.BCVdir;
                datasetV.Datasets_file      = app.Datasets_file;
                datasetV.Datasets           = app.Datasets;
                tableData.Name              = {app.Datasets.Dataset_name}';
                tableData.Description       = {app.Datasets.Description}';
                tableData.Path_Location     = {app.Datasets.Path}';
                datasetV.UITable.Data       = struct2table(tableData);
            end
        end

        % Callback function: ListDatasetPushTool
        function ListDatasetPushToolClicked(app, event)
            DatasetListMenuSelected(app, event);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = fileparts(mfilename('fullpath'));

            % Create BCVARETAToolboxv10UIFigure and hide until all components are created
            app.BCVARETAToolboxv10UIFigure = uifigure('Visible', 'off');
            app.BCVARETAToolboxv10UIFigure.Color = [0.9412 0.9412 0.9412];
            colormap(app.BCVARETAToolboxv10UIFigure, 'parula');
            app.BCVARETAToolboxv10UIFigure.Position = [100 100 698 517];
            app.BCVARETAToolboxv10UIFigure.Name = 'BC-VARETA Toolbox v1.0';
            app.BCVARETAToolboxv10UIFigure.WindowStyle = 'alwaysontop';

            % Create FileMenu
            app.FileMenu = uimenu(app.BCVARETAToolboxv10UIFigure);
            app.FileMenu.Text = 'File';

            % Create NewprocessMenu
            app.NewprocessMenu = uimenu(app.FileMenu);
            app.NewprocessMenu.Text = 'New process';

            % Create TestdataMenu
            app.TestdataMenu = uimenu(app.FileMenu);
            app.TestdataMenu.Text = 'Test data';

            % Create TestDataDownloadMenu
            app.TestDataDownloadMenu = uimenu(app.TestdataMenu);
            app.TestDataDownloadMenu.MenuSelectedFcn = createCallbackFcn(app, @TestDataDownloadMenuSelected, true);
            app.TestDataDownloadMenu.Text = 'Download';

            % Create TestDataImportMenu
            app.TestDataImportMenu = uimenu(app.TestdataMenu);
            app.TestDataImportMenu.MenuSelectedFcn = createCallbackFcn(app, @TestDataImportMenuSelected, true);
            app.TestDataImportMenu.Text = 'Import';

            % Create TestDataRunMenu
            app.TestDataRunMenu = uimenu(app.TestdataMenu);
            app.TestDataRunMenu.Text = 'Run';

            % Create DatasetMenu
            app.DatasetMenu = uimenu(app.FileMenu);
            app.DatasetMenu.Text = 'Dataset';

            % Create DatasetLoadMenu
            app.DatasetLoadMenu = uimenu(app.DatasetMenu);
            app.DatasetLoadMenu.MenuSelectedFcn = createCallbackFcn(app, @DatasetLoadMenuSelected, true);
            app.DatasetLoadMenu.Text = 'Load';

            % Create DatasetListMenu
            app.DatasetListMenu = uimenu(app.DatasetMenu);
            app.DatasetListMenu.MenuSelectedFcn = createCallbackFcn(app, @DatasetListMenuSelected, true);
            app.DatasetListMenu.Text = 'List';

            % Create DatasetDeleteMenu
            app.DatasetDeleteMenu = uimenu(app.DatasetMenu);
            app.DatasetDeleteMenu.Text = 'Delete';

            % Create ExitMenu
            app.ExitMenu = uimenu(app.FileMenu);
            app.ExitMenu.MenuSelectedFcn = createCallbackFcn(app, @ExitMenuSelected, true);
            app.ExitMenu.Text = 'Exit';

            % Create ToolsMenu
            app.ToolsMenu = uimenu(app.BCVARETAToolboxv10UIFigure);
            app.ToolsMenu.Text = 'Tools';

            % Create GroupstastMenu
            app.GroupstastMenu = uimenu(app.ToolsMenu);
            app.GroupstastMenu.Text = 'Group stast';

            % Create ViewMenu
            app.ViewMenu = uimenu(app.BCVARETAToolboxv10UIFigure);
            app.ViewMenu.Text = 'View';

            % Create DatasetresultsMenu
            app.DatasetresultsMenu = uimenu(app.ViewMenu);
            app.DatasetresultsMenu.Text = 'Dataset results';

            % Create HelpMenu
            app.HelpMenu = uimenu(app.BCVARETAToolboxv10UIFigure);
            app.HelpMenu.Text = 'Help';

            % Create CheckforupdateMenu
            app.CheckforupdateMenu = uimenu(app.HelpMenu);
            app.CheckforupdateMenu.Text = 'Check for update';

            % Create DocumentationMenu
            app.DocumentationMenu = uimenu(app.HelpMenu);
            app.DocumentationMenu.Text = 'Documentation';

            % Create AboutMenu
            app.AboutMenu = uimenu(app.HelpMenu);
            app.AboutMenu.Text = 'About';

            % Create Toolbar
            app.Toolbar = uitoolbar(app.BCVARETAToolboxv10UIFigure);

            % Create NewProcessPushTool
            app.NewProcessPushTool = uipushtool(app.Toolbar);
            app.NewProcessPushTool.Tooltip = {'New process'};
            app.NewProcessPushTool.Icon = fullfile(pathToMLAPP, 'images', 'new_process.jpeg');

            % Create Default_ParamsPushTool
            app.Default_ParamsPushTool = uipushtool(app.Toolbar);
            app.Default_ParamsPushTool.Tooltip = {'Set default params'};
            app.Default_ParamsPushTool.ClickedCallback = createCallbackFcn(app, @Default_ParamsPushToolClicked, true);
            app.Default_ParamsPushTool.Icon = fullfile(pathToMLAPP, 'images', 'default.png');

            % Create ListDatasetPushTool
            app.ListDatasetPushTool = uipushtool(app.Toolbar);
            app.ListDatasetPushTool.Tooltip = {'List Datasets'};
            app.ListDatasetPushTool.ClickedCallback = createCallbackFcn(app, @ListDatasetPushToolClicked, true);
            app.ListDatasetPushTool.Icon = fullfile(pathToMLAPP, 'images', 'list_dataset.png');

            % Create TabGroup
            app.TabGroup = uitabgroup(app.BCVARETAToolboxv10UIFigure);
            app.TabGroup.Position = [1 1 698 514];

            % Create PrincipalTab
            app.PrincipalTab = uitab(app.TabGroup);
            app.PrincipalTab.Title = 'Principal';

            % Create GeneralTab
            app.GeneralTab = uitab(app.TabGroup);
            app.GeneralTab.Tooltip = {'General parameters for the analysis'};
            app.GeneralTab.Title = 'General';

            % Create gener_PrincipalParamsPanel
            app.gener_PrincipalParamsPanel = uipanel(app.GeneralTab);
            app.gener_PrincipalParamsPanel.Title = 'Principal params';
            app.gener_PrincipalParamsPanel.FontWeight = 'bold';
            app.gener_PrincipalParamsPanel.FontSize = 14;
            app.gener_PrincipalParamsPanel.Position = [11 52 675 426];

            % Create gener_SelectOutputButton
            app.gener_SelectOutputButton = uibutton(app.gener_PrincipalParamsPanel, 'push');
            app.gener_SelectOutputButton.ButtonPushedFcn = createCallbackFcn(app, @gener_SelectOutputButtonPushed, true);
            app.gener_SelectOutputButton.FontSize = 14;
            app.gener_SelectOutputButton.Position = [591 147 68 25];
            app.gener_SelectOutputButton.Text = 'Select';

            % Create gener_SelectInputButton
            app.gener_SelectInputButton = uibutton(app.gener_PrincipalParamsPanel, 'push');
            app.gener_SelectInputButton.ButtonPushedFcn = createCallbackFcn(app, @gener_SelectInputButtonPushed, true);
            app.gener_SelectInputButton.FontSize = 14;
            app.gener_SelectInputButton.Position = [591 182 68 25];
            app.gener_SelectInputButton.Text = 'Select';

            % Create gener_InputfolderEditField
            app.gener_InputfolderEditField = uieditfield(app.gener_PrincipalParamsPanel, 'text');
            app.gener_InputfolderEditField.FontSize = 14;
            app.gener_InputfolderEditField.Tooltip = {'Input path adderess wich include all the subjects to analysis.'};
            app.gener_InputfolderEditField.Position = [112 185 473 22];

            % Create InputfolderLabel
            app.InputfolderLabel = uilabel(app.gener_PrincipalParamsPanel);
            app.InputfolderLabel.HorizontalAlignment = 'right';
            app.InputfolderLabel.FontSize = 14;
            app.InputfolderLabel.Position = [19 185 79 22];
            app.InputfolderLabel.Text = 'Input folder:';

            % Create gener_OutputfolderEditField
            app.gener_OutputfolderEditField = uieditfield(app.gener_PrincipalParamsPanel, 'text');
            app.gener_OutputfolderEditField.ValueChangedFcn = createCallbackFcn(app, @gener_OutputfolderEditFieldValueChanged, true);
            app.gener_OutputfolderEditField.FontSize = 14;
            app.gener_OutputfolderEditField.Tooltip = {'Output path adderess where will be saved the BC-VARETA results.'};
            app.gener_OutputfolderEditField.Position = [111 148 473 22];

            % Create OutputfolderLabel
            app.OutputfolderLabel = uilabel(app.gener_PrincipalParamsPanel);
            app.OutputfolderLabel.HorizontalAlignment = 'right';
            app.OutputfolderLabel.FontSize = 14;
            app.OutputfolderLabel.Position = [8 148 90 22];
            app.OutputfolderLabel.Text = 'Output folder:';

            % Create gener_AnalysislevelDropDown
            app.gener_AnalysislevelDropDown = uidropdown(app.gener_PrincipalParamsPanel);
            app.gener_AnalysislevelDropDown.Items = {'All analysis', 'Sensor', 'Activation', 'Connectivity', 'Sensor and Activation', 'Activation and Connectivity'};
            app.gener_AnalysislevelDropDown.ItemsData = {'all', '1', '2', '3', '12', '23'};
            app.gener_AnalysislevelDropDown.ValueChangedFcn = createCallbackFcn(app, @gener_AnalysislevelDropDownValueChanged, true);
            app.gener_AnalysislevelDropDown.Tooltip = {'1- All analysis (sensor, activation and connectivity)'; '2- Compute only sensor.'; '3- Compute only activation.'; '4- Compute only connectivity.'; '5- Compute sensor and activation.'; '6- Compute activation and connectivity.'};
            app.gener_AnalysislevelDropDown.FontSize = 14;
            app.gener_AnalysislevelDropDown.Position = [109 109 229 22];
            app.gener_AnalysislevelDropDown.Value = 'all';

            % Create AnalysislevelDropDownLabel
            app.AnalysislevelDropDownLabel = uilabel(app.gener_PrincipalParamsPanel);
            app.AnalysislevelDropDownLabel.HorizontalAlignment = 'right';
            app.AnalysislevelDropDownLabel.FontSize = 14;
            app.AnalysislevelDropDownLabel.Position = [4 109 94 22];
            app.AnalysislevelDropDownLabel.Text = 'Analysis level:';

            % Create gener_RunbytrialsSwitch
            app.gener_RunbytrialsSwitch = uiswitch(app.gener_PrincipalParamsPanel, 'slider');
            app.gener_RunbytrialsSwitch.Items = {'No', 'Yes'};
            app.gener_RunbytrialsSwitch.ValueChangedFcn = createCallbackFcn(app, @gener_RunbytrialsSwitchValueChanged, true);
            app.gener_RunbytrialsSwitch.Tooltip = {'Yes - The time series of different trials are analyzed independently'; 'No - The time series are combined in a single time serie.'};
            app.gener_RunbytrialsSwitch.FontSize = 14;
            app.gener_RunbytrialsSwitch.Position = [136 70 45 20];
            app.gener_RunbytrialsSwitch.Value = 'No';

            % Create RunbytrialsSwitchLabel
            app.RunbytrialsSwitchLabel = uilabel(app.gener_PrincipalParamsPanel);
            app.RunbytrialsSwitchLabel.HorizontalAlignment = 'center';
            app.RunbytrialsSwitchLabel.FontSize = 14;
            app.RunbytrialsSwitchLabel.Position = [13 70 87 22];
            app.RunbytrialsSwitchLabel.Text = 'Run by trials:';

            % Create gener_RunbyfrequencybinSwitch
            app.gener_RunbyfrequencybinSwitch = uiswitch(app.gener_PrincipalParamsPanel, 'slider');
            app.gener_RunbyfrequencybinSwitch.Items = {'No', 'Yes'};
            app.gener_RunbyfrequencybinSwitch.ValueChangedFcn = createCallbackFcn(app, @gener_RunbyfrequencybinSwitchValueChanged, true);
            app.gener_RunbyfrequencybinSwitch.Tooltip = {'Yes - Frequency components are analyzed independentlly'; 'No - According to a band structured sparsity model corresponding to smoothness in the spectral domain.'};
            app.gener_RunbyfrequencybinSwitch.FontSize = 14;
            app.gener_RunbyfrequencybinSwitch.Position = [437 68 45 20];
            app.gener_RunbyfrequencybinSwitch.Value = 'No';

            % Create RunbyfrequencybinSwitchLabel
            app.RunbyfrequencybinSwitchLabel = uilabel(app.gener_PrincipalParamsPanel);
            app.RunbyfrequencybinSwitchLabel.HorizontalAlignment = 'center';
            app.RunbyfrequencybinSwitchLabel.FontSize = 14;
            app.RunbyfrequencybinSwitchLabel.Position = [265 67 141 22];
            app.RunbyfrequencybinSwitchLabel.Text = 'Run by frequency bin:';

            % Create gener_GetsystemresponseSwitch
            app.gener_GetsystemresponseSwitch = uiswitch(app.gener_PrincipalParamsPanel, 'slider');
            app.gener_GetsystemresponseSwitch.Items = {'No', 'Yes'};
            app.gener_GetsystemresponseSwitch.ValueChangedFcn = createCallbackFcn(app, @gener_GetsystemresponseSwitchValueChanged, true);
            app.gener_GetsystemresponseSwitch.Tooltip = {'Decide wheather the souces are responsive for all frequencies or not when applying sSSBL z-scores.'; 'No - Apply z-scores indepoendentlly for all frequencies.'; 'Yes - Apply z-scores jointly for all frequencies.'};
            app.gener_GetsystemresponseSwitch.FontSize = 14;
            app.gener_GetsystemresponseSwitch.Position = [438 25 45 20];
            app.gener_GetsystemresponseSwitch.Value = 'No';

            % Create GetsystemresponseSwitchLabel
            app.GetsystemresponseSwitchLabel = uilabel(app.gener_PrincipalParamsPanel);
            app.GetsystemresponseSwitchLabel.HorizontalAlignment = 'center';
            app.GetsystemresponseSwitchLabel.FontSize = 14;
            app.GetsystemresponseSwitchLabel.Position = [266 24 141 22];
            app.GetsystemresponseSwitchLabel.Text = 'Get system response:';

            % Create RuninGPUSwitchLabel
            app.RuninGPUSwitchLabel = uilabel(app.gener_PrincipalParamsPanel);
            app.RuninGPUSwitchLabel.HorizontalAlignment = 'center';
            app.RuninGPUSwitchLabel.FontSize = 14;
            app.RuninGPUSwitchLabel.Position = [16 25 84 22];
            app.RuninGPUSwitchLabel.Text = 'Run in GPU:';

            % Create gener_RuninGPUSwitch
            app.gener_RuninGPUSwitch = uiswitch(app.gener_PrincipalParamsPanel, 'slider');
            app.gener_RuninGPUSwitch.Items = {'No', 'Yes'};
            app.gener_RuninGPUSwitch.ValueChangedFcn = createCallbackFcn(app, @gener_RuninGPUSwitchValueChanged, true);
            app.gener_RuninGPUSwitch.Tooltip = {'Yes - Use GPU in the processing.'; 'No - Do not use GPU in the processing.'};
            app.gener_RuninGPUSwitch.FontSize = 14;
            app.gener_RuninGPUSwitch.Position = [137 26 45 20];
            app.gener_RuninGPUSwitch.Value = 'No';

            % Create gener_DatasetNameEditField
            app.gener_DatasetNameEditField = uieditfield(app.gener_PrincipalParamsPanel, 'text');
            app.gener_DatasetNameEditField.FontSize = 14;
            app.gener_DatasetNameEditField.Tooltip = {'Input path adderess wich include all the subjects to analysis.'};
            app.gener_DatasetNameEditField.Position = [113 363 217 22];

            % Create InputfolderLabel_2
            app.InputfolderLabel_2 = uilabel(app.gener_PrincipalParamsPanel);
            app.InputfolderLabel_2.HorizontalAlignment = 'right';
            app.InputfolderLabel_2.FontSize = 14;
            app.InputfolderLabel_2.Position = [3 363 96 22];
            app.InputfolderLabel_2.Text = 'Dataset name:';

            % Create DescriptionTextAreaLabel
            app.DescriptionTextAreaLabel = uilabel(app.gener_PrincipalParamsPanel);
            app.DescriptionTextAreaLabel.HorizontalAlignment = 'right';
            app.DescriptionTextAreaLabel.FontSize = 14;
            app.DescriptionTextAreaLabel.Position = [19 324 79 22];
            app.DescriptionTextAreaLabel.Text = 'Description:';

            % Create gener_DescriptionTextArea
            app.gener_DescriptionTextArea = uitextarea(app.gener_PrincipalParamsPanel);
            app.gener_DescriptionTextArea.FontSize = 14;
            app.gener_DescriptionTextArea.Position = [113 220 546 128];

            % Create SensorTab
            app.SensorTab = uitab(app.TabGroup);
            app.SensorTab.Tooltip = {'Parameters in the first level of analysis or sensor.'};
            app.SensorTab.Title = 'Sensor';

            % Create spect_FreqPanel
            app.spect_FreqPanel = uipanel(app.SensorTab);
            app.spect_FreqPanel.Tooltip = {'Definition of frequency bands for analysis.'};
            app.spect_FreqPanel.Title = 'Frequency bands';
            app.spect_FreqPanel.FontWeight = 'bold';
            app.spect_FreqPanel.FontSize = 14;
            app.spect_FreqPanel.Position = [6 51 680 258];

            % Create GridLayout
            app.GridLayout = uigridlayout(app.spect_FreqPanel);
            app.GridLayout.ColumnWidth = {'1x', '1x', '1x', '1x'};
            app.GridLayout.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout.ColumnSpacing = 5;
            app.GridLayout.RowSpacing = 2;
            app.GridLayout.Padding = [5 5 5 5];

            % Create NameLabel
            app.NameLabel = uilabel(app.GridLayout);
            app.NameLabel.HorizontalAlignment = 'center';
            app.NameLabel.FontSize = 14;
            app.NameLabel.FontWeight = 'bold';
            app.NameLabel.Layout.Row = 1;
            app.NameLabel.Layout.Column = 1;
            app.NameLabel.Text = 'Name';

            % Create ProcessLabel
            app.ProcessLabel = uilabel(app.GridLayout);
            app.ProcessLabel.HorizontalAlignment = 'center';
            app.ProcessLabel.FontSize = 14;
            app.ProcessLabel.FontWeight = 'bold';
            app.ProcessLabel.Layout.Row = 1;
            app.ProcessLabel.Layout.Column = 2;
            app.ProcessLabel.Text = 'Process';

            % Create StartLabel
            app.StartLabel = uilabel(app.GridLayout);
            app.StartLabel.HorizontalAlignment = 'center';
            app.StartLabel.FontSize = 14;
            app.StartLabel.FontWeight = 'bold';
            app.StartLabel.Layout.Row = 1;
            app.StartLabel.Layout.Column = 3;
            app.StartLabel.Text = 'Start';

            % Create EndLabel
            app.EndLabel = uilabel(app.GridLayout);
            app.EndLabel.HorizontalAlignment = 'center';
            app.EndLabel.FontSize = 14;
            app.EndLabel.FontWeight = 'bold';
            app.EndLabel.Layout.Row = 1;
            app.EndLabel.Layout.Column = 4;
            app.EndLabel.Text = 'End';

            % Create delta_startEditField
            app.delta_startEditField = uieditfield(app.GridLayout, 'numeric');
            app.delta_startEditField.Limits = [0 4];
            app.delta_startEditField.ValueChangedFcn = createCallbackFcn(app, @delta_startEditFieldValueChanged, true);
            app.delta_startEditField.FontSize = 14;
            app.delta_startEditField.Layout.Row = 2;
            app.delta_startEditField.Layout.Column = 3;

            % Create delta_endEditField
            app.delta_endEditField = uieditfield(app.GridLayout, 'numeric');
            app.delta_endEditField.Limits = [0 4];
            app.delta_endEditField.ValueChangedFcn = createCallbackFcn(app, @delta_endEditFieldValueChanged, true);
            app.delta_endEditField.FontSize = 14;
            app.delta_endEditField.Layout.Row = 2;
            app.delta_endEditField.Layout.Column = 4;
            app.delta_endEditField.Value = 4;

            % Create theta_startEditField
            app.theta_startEditField = uieditfield(app.GridLayout, 'numeric');
            app.theta_startEditField.Limits = [4 7];
            app.theta_startEditField.ValueChangedFcn = createCallbackFcn(app, @theta_startEditFieldValueChanged, true);
            app.theta_startEditField.FontSize = 14;
            app.theta_startEditField.Layout.Row = 3;
            app.theta_startEditField.Layout.Column = 3;
            app.theta_startEditField.Value = 4;

            % Create theta_endEditField
            app.theta_endEditField = uieditfield(app.GridLayout, 'numeric');
            app.theta_endEditField.Limits = [4 7];
            app.theta_endEditField.ValueChangedFcn = createCallbackFcn(app, @theta_endEditFieldValueChanged, true);
            app.theta_endEditField.FontSize = 14;
            app.theta_endEditField.Layout.Row = 3;
            app.theta_endEditField.Layout.Column = 4;
            app.theta_endEditField.Value = 7;

            % Create alpha_startEditField
            app.alpha_startEditField = uieditfield(app.GridLayout, 'numeric');
            app.alpha_startEditField.Limits = [7 14];
            app.alpha_startEditField.ValueChangedFcn = createCallbackFcn(app, @alpha_startEditFieldValueChanged, true);
            app.alpha_startEditField.FontSize = 14;
            app.alpha_startEditField.Layout.Row = 4;
            app.alpha_startEditField.Layout.Column = 3;
            app.alpha_startEditField.Value = 7;

            % Create alpha_endEditField
            app.alpha_endEditField = uieditfield(app.GridLayout, 'numeric');
            app.alpha_endEditField.Limits = [7 14];
            app.alpha_endEditField.ValueChangedFcn = createCallbackFcn(app, @alpha_endEditFieldValueChanged, true);
            app.alpha_endEditField.FontSize = 14;
            app.alpha_endEditField.Layout.Row = 4;
            app.alpha_endEditField.Layout.Column = 4;
            app.alpha_endEditField.Value = 14;

            % Create beta_startEditField
            app.beta_startEditField = uieditfield(app.GridLayout, 'numeric');
            app.beta_startEditField.Limits = [14 31];
            app.beta_startEditField.ValueChangedFcn = createCallbackFcn(app, @beta_startEditFieldValueChanged, true);
            app.beta_startEditField.FontSize = 14;
            app.beta_startEditField.Layout.Row = 5;
            app.beta_startEditField.Layout.Column = 3;
            app.beta_startEditField.Value = 14;

            % Create beta_endEditField
            app.beta_endEditField = uieditfield(app.GridLayout, 'numeric');
            app.beta_endEditField.Limits = [14 31];
            app.beta_endEditField.ValueChangedFcn = createCallbackFcn(app, @beta_endEditFieldValueChanged, true);
            app.beta_endEditField.FontSize = 14;
            app.beta_endEditField.Layout.Row = 5;
            app.beta_endEditField.Layout.Column = 4;
            app.beta_endEditField.Value = 31;

            % Create gamma1_startEditField
            app.gamma1_startEditField = uieditfield(app.GridLayout, 'numeric');
            app.gamma1_startEditField.Limits = [31 60];
            app.gamma1_startEditField.ValueChangedFcn = createCallbackFcn(app, @gamma1_startEditFieldValueChanged, true);
            app.gamma1_startEditField.FontSize = 14;
            app.gamma1_startEditField.Layout.Row = 6;
            app.gamma1_startEditField.Layout.Column = 3;
            app.gamma1_startEditField.Value = 31;

            % Create gamma1_endEditField
            app.gamma1_endEditField = uieditfield(app.GridLayout, 'numeric');
            app.gamma1_endEditField.Limits = [31 60];
            app.gamma1_endEditField.ValueChangedFcn = createCallbackFcn(app, @gamma1_endEditFieldValueChanged, true);
            app.gamma1_endEditField.FontSize = 14;
            app.gamma1_endEditField.Layout.Row = 6;
            app.gamma1_endEditField.Layout.Column = 4;
            app.gamma1_endEditField.Value = 60;

            % Create gamma2_startEditField
            app.gamma2_startEditField = uieditfield(app.GridLayout, 'numeric');
            app.gamma2_startEditField.Limits = [60 90];
            app.gamma2_startEditField.ValueChangedFcn = createCallbackFcn(app, @gamma2_startEditFieldValueChanged, true);
            app.gamma2_startEditField.FontSize = 14;
            app.gamma2_startEditField.Layout.Row = 7;
            app.gamma2_startEditField.Layout.Column = 3;
            app.gamma2_startEditField.Value = 60;

            % Create gamma2_endEditField
            app.gamma2_endEditField = uieditfield(app.GridLayout, 'numeric');
            app.gamma2_endEditField.Limits = [60 90];
            app.gamma2_endEditField.ValueChangedFcn = createCallbackFcn(app, @gamma2_endEditFieldValueChanged, true);
            app.gamma2_endEditField.FontSize = 14;
            app.gamma2_endEditField.Layout.Row = 7;
            app.gamma2_endEditField.Layout.Column = 4;
            app.gamma2_endEditField.Value = 90;

            % Create delta_runSwitch
            app.delta_runSwitch = uiswitch(app.GridLayout, 'slider');
            app.delta_runSwitch.Items = {'No', 'Yes'};
            app.delta_runSwitch.ValueChangedFcn = createCallbackFcn(app, @delta_runSwitchValueChanged, true);
            app.delta_runSwitch.FontSize = 14;
            app.delta_runSwitch.Layout.Row = 2;
            app.delta_runSwitch.Layout.Column = 2;
            app.delta_runSwitch.Value = 'Yes';

            % Create theta_runSwitch
            app.theta_runSwitch = uiswitch(app.GridLayout, 'slider');
            app.theta_runSwitch.Items = {'No', 'Yes'};
            app.theta_runSwitch.ValueChangedFcn = createCallbackFcn(app, @theta_runSwitchValueChanged, true);
            app.theta_runSwitch.FontSize = 14;
            app.theta_runSwitch.Layout.Row = 3;
            app.theta_runSwitch.Layout.Column = 2;
            app.theta_runSwitch.Value = 'Yes';

            % Create alpha_runSwitch
            app.alpha_runSwitch = uiswitch(app.GridLayout, 'slider');
            app.alpha_runSwitch.Items = {'No', 'Yes'};
            app.alpha_runSwitch.ValueChangedFcn = createCallbackFcn(app, @alpha_runSwitchValueChanged, true);
            app.alpha_runSwitch.FontSize = 14;
            app.alpha_runSwitch.Layout.Row = 4;
            app.alpha_runSwitch.Layout.Column = 2;
            app.alpha_runSwitch.Value = 'Yes';

            % Create beta_runSwitch
            app.beta_runSwitch = uiswitch(app.GridLayout, 'slider');
            app.beta_runSwitch.Items = {'No', 'Yes'};
            app.beta_runSwitch.ValueChangedFcn = createCallbackFcn(app, @beta_runSwitchValueChanged, true);
            app.beta_runSwitch.FontSize = 14;
            app.beta_runSwitch.Layout.Row = 5;
            app.beta_runSwitch.Layout.Column = 2;
            app.beta_runSwitch.Value = 'Yes';

            % Create gamma1_runSwitch
            app.gamma1_runSwitch = uiswitch(app.GridLayout, 'slider');
            app.gamma1_runSwitch.Items = {'No', 'Yes'};
            app.gamma1_runSwitch.ValueChangedFcn = createCallbackFcn(app, @gamma1_runSwitchValueChanged, true);
            app.gamma1_runSwitch.FontSize = 14;
            app.gamma1_runSwitch.Layout.Row = 6;
            app.gamma1_runSwitch.Layout.Column = 2;
            app.gamma1_runSwitch.Value = 'Yes';

            % Create gamma2_runSwitch
            app.gamma2_runSwitch = uiswitch(app.GridLayout, 'slider');
            app.gamma2_runSwitch.Items = {'No', 'Yes'};
            app.gamma2_runSwitch.ValueChangedFcn = createCallbackFcn(app, @gamma2_runSwitchValueChanged, true);
            app.gamma2_runSwitch.FontSize = 14;
            app.gamma2_runSwitch.Layout.Row = 7;
            app.gamma2_runSwitch.Layout.Column = 2;
            app.gamma2_runSwitch.Value = 'Yes';

            % Create DeltaLabel
            app.DeltaLabel = uilabel(app.GridLayout);
            app.DeltaLabel.HorizontalAlignment = 'center';
            app.DeltaLabel.FontSize = 14;
            app.DeltaLabel.Tooltip = {'delta band. the starting and ending frequencies of the band are start and end. set run in:'; 'Yes -  if you want this band to be analyzed'};
            app.DeltaLabel.Layout.Row = 2;
            app.DeltaLabel.Layout.Column = 1;
            app.DeltaLabel.Text = 'Delta';

            % Create ThetaLabel
            app.ThetaLabel = uilabel(app.GridLayout);
            app.ThetaLabel.HorizontalAlignment = 'center';
            app.ThetaLabel.FontSize = 14;
            app.ThetaLabel.Tooltip = {'theta band. the starting and ending frequencies of the band are start and end. set run in:'; 'Yes -  if you want this band to be analyzed'};
            app.ThetaLabel.Layout.Row = 3;
            app.ThetaLabel.Layout.Column = 1;
            app.ThetaLabel.Text = 'Theta';

            % Create AlphaLabel
            app.AlphaLabel = uilabel(app.GridLayout);
            app.AlphaLabel.HorizontalAlignment = 'center';
            app.AlphaLabel.FontSize = 14;
            app.AlphaLabel.Tooltip = {'alpha band. the starting and ending frequencies of the band are start and end. set run in:'; 'Yes -  if you want this band to be analyzed'};
            app.AlphaLabel.Layout.Row = 4;
            app.AlphaLabel.Layout.Column = 1;
            app.AlphaLabel.Text = 'Alpha';

            % Create BetaLabel
            app.BetaLabel = uilabel(app.GridLayout);
            app.BetaLabel.HorizontalAlignment = 'center';
            app.BetaLabel.FontSize = 14;
            app.BetaLabel.Tooltip = {'beta band. the starting and ending frequencies of the band are start and end. set run in:'; 'Yes -  if you want this band to be analyzed'};
            app.BetaLabel.Layout.Row = 5;
            app.BetaLabel.Layout.Column = 1;
            app.BetaLabel.Text = 'Beta';

            % Create Gamma1Label
            app.Gamma1Label = uilabel(app.GridLayout);
            app.Gamma1Label.HorizontalAlignment = 'center';
            app.Gamma1Label.FontSize = 14;
            app.Gamma1Label.Tooltip = {'gamma1 band. the starting and ending frequencies of the band are start and end. set run in:'; 'Yes -  if you want this band to be analyzed'};
            app.Gamma1Label.Layout.Row = 6;
            app.Gamma1Label.Layout.Column = 1;
            app.Gamma1Label.Text = 'Gamma1';

            % Create Gamma2Label
            app.Gamma2Label = uilabel(app.GridLayout);
            app.Gamma2Label.HorizontalAlignment = 'center';
            app.Gamma2Label.FontSize = 14;
            app.Gamma2Label.Tooltip = {'gamma2 band. the starting and ending frequencies of the band are start and end. set run in:'; 'Yes -  if you want this band to be analyzed'};
            app.Gamma2Label.Layout.Row = 7;
            app.Gamma2Label.Layout.Column = 1;
            app.Gamma2Label.Text = 'Gamma2';

            % Create spect_OtherparamsPanel
            app.spect_OtherparamsPanel = uipanel(app.SensorTab);
            app.spect_OtherparamsPanel.Title = 'Data params';
            app.spect_OtherparamsPanel.FontWeight = 'bold';
            app.spect_OtherparamsPanel.FontSize = 14;
            app.spect_OtherparamsPanel.Position = [6 316 680 170];

            % Create FreqresolutionLabel
            app.FreqresolutionLabel = uilabel(app.spect_OtherparamsPanel);
            app.FreqresolutionLabel.HorizontalAlignment = 'right';
            app.FreqresolutionLabel.FontSize = 14;
            app.FreqresolutionLabel.Position = [11 62 102 22];
            app.FreqresolutionLabel.Text = 'Freq resolution:';

            % Create spect_freq_resolEditField
            app.spect_freq_resolEditField = uieditfield(app.spect_OtherparamsPanel, 'numeric');
            app.spect_freq_resolEditField.Limits = [0.1 1];
            app.spect_freq_resolEditField.FontSize = 14;
            app.spect_freq_resolEditField.Tooltip = {'Spectral resolution of the FFT'};
            app.spect_freq_resolEditField.Position = [121 62 67 22];
            app.spect_freq_resolEditField.Value = 0.5;

            % Create SamplefreqLabel
            app.SamplefreqLabel = uilabel(app.spect_OtherparamsPanel);
            app.SamplefreqLabel.HorizontalAlignment = 'right';
            app.SamplefreqLabel.FontSize = 14;
            app.SamplefreqLabel.Position = [256 62 84 22];
            app.SamplefreqLabel.Text = 'Sample freq:';

            % Create spect_sample_freqEditField
            app.spect_sample_freqEditField = uieditfield(app.spect_OtherparamsPanel, 'numeric');
            app.spect_sample_freqEditField.Limits = [100 600];
            app.spect_sample_freqEditField.FontSize = 14;
            app.spect_sample_freqEditField.Tooltip = {'MEG/EEG sampling frequency'};
            app.spect_sample_freqEditField.Position = [346 62 59 22];
            app.spect_sample_freqEditField.Value = 200;

            % Create FreqmaxEditFieldLabel
            app.FreqmaxEditFieldLabel = uilabel(app.spect_OtherparamsPanel);
            app.FreqmaxEditFieldLabel.HorizontalAlignment = 'right';
            app.FreqmaxEditFieldLabel.FontSize = 14;
            app.FreqmaxEditFieldLabel.Position = [271 20 68 22];
            app.FreqmaxEditFieldLabel.Text = 'Freq max:';

            % Create spect_FreqmaxEditField
            app.spect_FreqmaxEditField = uieditfield(app.spect_OtherparamsPanel, 'numeric');
            app.spect_FreqmaxEditField.Limits = [0 90];
            app.spect_FreqmaxEditField.FontSize = 14;
            app.spect_FreqmaxEditField.Tooltip = {'Maximum frequency under analysis. equal or higher than the ending of the faster band considered'};
            app.spect_FreqmaxEditField.Position = [347 20 57 22];
            app.spect_FreqmaxEditField.Value = 90;

            % Create freq_gfiltvarEditFieldLabel
            app.freq_gfiltvarEditFieldLabel = uilabel(app.spect_OtherparamsPanel);
            app.freq_gfiltvarEditFieldLabel.HorizontalAlignment = 'right';
            app.freq_gfiltvarEditFieldLabel.FontSize = 14;
            app.freq_gfiltvarEditFieldLabel.Position = [31 21 82 22];
            app.freq_gfiltvarEditFieldLabel.Text = 'freq_gfiltvar:';

            % Create spect_freq_gfiltvarEditField
            app.spect_freq_gfiltvarEditField = uieditfield(app.spect_OtherparamsPanel, 'numeric');
            app.spect_freq_gfiltvarEditField.Limits = [0.1 3];
            app.spect_freq_gfiltvarEditField.FontSize = 14;
            app.spect_freq_gfiltvarEditField.Tooltip = {'Variance (given in Hz) for the gaussian response filter of the data'};
            app.spect_freq_gfiltvarEditField.Position = [122 21 69 22];
            app.spect_freq_gfiltvarEditField.Value = 0.5;

            % Create win_orderLabel
            app.win_orderLabel = uilabel(app.spect_OtherparamsPanel);
            app.win_orderLabel.HorizontalAlignment = 'right';
            app.win_orderLabel.FontSize = 14;
            app.win_orderLabel.Position = [468 63 70 22];
            app.win_orderLabel.Text = 'win_order:';

            % Create spect_win_orderSlider
            app.spect_win_orderSlider = uislider(app.spect_OtherparamsPanel);
            app.spect_win_orderSlider.Limits = [1 3];
            app.spect_win_orderSlider.MajorTicks = [1 2 3];
            app.spect_win_orderSlider.MinorTicks = [1 2 3];
            app.spect_win_orderSlider.Tooltip = {'Windowing Order. number of Slepian sequences in this case'};
            app.spect_win_orderSlider.FontSize = 14;
            app.spect_win_orderSlider.Position = [560 74 74 3];
            app.spect_win_orderSlider.Value = 1;

            % Create FTTmethodLabel
            app.FTTmethodLabel = uilabel(app.spect_OtherparamsPanel);
            app.FTTmethodLabel.HorizontalAlignment = 'right';
            app.FTTmethodLabel.FontSize = 14;
            app.FTTmethodLabel.Position = [29 106 85 22];
            app.FTTmethodLabel.Text = 'FFT method:';

            % Create sens_FFTMethodDropDown
            app.sens_FFTMethodDropDown = uidropdown(app.spect_OtherparamsPanel);
            app.sens_FFTMethodDropDown.Items = {'Hilbert', 'Thomson'};
            app.sens_FFTMethodDropDown.FontSize = 14;
            app.sens_FFTMethodDropDown.Position = [122 106 137 22];
            app.sens_FFTMethodDropDown.Value = 'Hilbert';

            % Create ActivationTab
            app.ActivationTab = uitab(app.TabGroup);
            app.ActivationTab.Tooltip = {'Parameters in the second level of analysis or activation.'};
            app.ActivationTab.Title = 'Activation';

            % Create activ_MethodsPanel
            app.activ_MethodsPanel = uipanel(app.ActivationTab);
            app.activ_MethodsPanel.Title = 'Methods';
            app.activ_MethodsPanel.FontWeight = 'bold';
            app.activ_MethodsPanel.FontSize = 14;
            app.activ_MethodsPanel.Position = [9 235 676 246];

            % Create SSSBLppLabel
            app.SSSBLppLabel = uilabel(app.activ_MethodsPanel);
            app.SSSBLppLabel.FontSize = 14;
            app.SSSBLppLabel.FontWeight = 'bold';
            app.SSSBLppLabel.Position = [23 193 73 22];
            app.SSSBLppLabel.Text = 'SSSBLpp:';

            % Create CheckBox_sssblpp
            app.CheckBox_sssblpp = uicheckbox(app.activ_MethodsPanel);
            app.CheckBox_sssblpp.Text = '';
            app.CheckBox_sssblpp.Position = [98 192 25 22];

            % Create eLORETALabel
            app.eLORETALabel = uilabel(app.activ_MethodsPanel);
            app.eLORETALabel.FontSize = 14;
            app.eLORETALabel.FontWeight = 'bold';
            app.eLORETALabel.Position = [235 193 74 22];
            app.eLORETALabel.Text = 'eLORETA:';

            % Create CheckBox_eloreta
            app.CheckBox_eloreta = uicheckbox(app.activ_MethodsPanel);
            app.CheckBox_eloreta.Text = '';
            app.CheckBox_eloreta.Position = [312 193 25 22];

            % Create LCMVLabel
            app.LCMVLabel = uilabel(app.activ_MethodsPanel);
            app.LCMVLabel.FontSize = 14;
            app.LCMVLabel.FontWeight = 'bold';
            app.LCMVLabel.Position = [473 193 49 22];
            app.LCMVLabel.Text = 'LCMV:';

            % Create CheckBox_lcmv
            app.CheckBox_lcmv = uicheckbox(app.activ_MethodsPanel);
            app.CheckBox_lcmv.Text = '';
            app.CheckBox_lcmv.Position = [529 193 25 22];

            % Create ThresholdEditFieldLabel
            app.ThresholdEditFieldLabel = uilabel(app.activ_MethodsPanel);
            app.ThresholdEditFieldLabel.HorizontalAlignment = 'right';
            app.ThresholdEditFieldLabel.FontSize = 14;
            app.ThresholdEditFieldLabel.Position = [19 155 71 22];
            app.ThresholdEditFieldLabel.Text = 'Threshold:';

            % Create sssblpp_thEditField
            app.sssblpp_thEditField = uieditfield(app.activ_MethodsPanel, 'numeric');
            app.sssblpp_thEditField.Limits = [0 1.7321];
            app.sssblpp_thEditField.FontSize = 14;
            app.sssblpp_thEditField.Position = [99 155 54 22];
            app.sssblpp_thEditField.Value = 1;

            % Create ThresholdEditField_2Label
            app.ThresholdEditField_2Label = uilabel(app.activ_MethodsPanel);
            app.ThresholdEditField_2Label.HorizontalAlignment = 'right';
            app.ThresholdEditField_2Label.FontSize = 14;
            app.ThresholdEditField_2Label.Position = [231 155 71 22];
            app.ThresholdEditField_2Label.Text = 'Threshold:';

            % Create eloreta_thEditField
            app.eloreta_thEditField = uieditfield(app.activ_MethodsPanel, 'numeric');
            app.eloreta_thEditField.Limits = [0 1.7321];
            app.eloreta_thEditField.FontSize = 14;
            app.eloreta_thEditField.Position = [312 155 53 22];
            app.eloreta_thEditField.Value = 1;

            % Create gamma1EditFieldLabel
            app.gamma1EditFieldLabel = uilabel(app.activ_MethodsPanel);
            app.gamma1EditFieldLabel.HorizontalAlignment = 'right';
            app.gamma1EditFieldLabel.FontSize = 14;
            app.gamma1EditFieldLabel.Position = [239 115 63 22];
            app.gamma1EditFieldLabel.Text = 'gamma1:';

            % Create eloreta_gamma1EditField
            app.eloreta_gamma1EditField = uieditfield(app.activ_MethodsPanel, 'numeric');
            app.eloreta_gamma1EditField.Limits = [0 2];
            app.eloreta_gamma1EditField.FontSize = 14;
            app.eloreta_gamma1EditField.Position = [312 115 53 22];

            % Create gamma2EditFieldLabel
            app.gamma2EditFieldLabel = uilabel(app.activ_MethodsPanel);
            app.gamma2EditFieldLabel.HorizontalAlignment = 'right';
            app.gamma2EditFieldLabel.FontSize = 14;
            app.gamma2EditFieldLabel.Position = [239 76 63 22];
            app.gamma2EditFieldLabel.Text = 'gamma2:';

            % Create eloreta_gamma2EditField
            app.eloreta_gamma2EditField = uieditfield(app.activ_MethodsPanel, 'numeric');
            app.eloreta_gamma2EditField.Limits = [0 2];
            app.eloreta_gamma2EditField.FontSize = 14;
            app.eloreta_gamma2EditField.Position = [313 76 53 22];
            app.eloreta_gamma2EditField.Value = 2;

            % Create deltagammaEditFieldLabel
            app.deltagammaEditFieldLabel = uilabel(app.activ_MethodsPanel);
            app.deltagammaEditFieldLabel.HorizontalAlignment = 'right';
            app.deltagammaEditFieldLabel.FontSize = 14;
            app.deltagammaEditFieldLabel.Position = [212 37 90 22];
            app.deltagammaEditFieldLabel.Text = 'delta gamma:';

            % Create eloreta_deltagammaEditField
            app.eloreta_deltagammaEditField = uieditfield(app.activ_MethodsPanel, 'numeric');
            app.eloreta_deltagammaEditField.Limits = [0.1 1];
            app.eloreta_deltagammaEditField.FontSize = 14;
            app.eloreta_deltagammaEditField.Position = [314 37 53 22];
            app.eloreta_deltagammaEditField.Value = 0.1;

            % Create ThresholdEditField_3Label
            app.ThresholdEditField_3Label = uilabel(app.activ_MethodsPanel);
            app.ThresholdEditField_3Label.HorizontalAlignment = 'right';
            app.ThresholdEditField_3Label.FontSize = 14;
            app.ThresholdEditField_3Label.Position = [447 155 71 22];
            app.ThresholdEditField_3Label.Text = 'Threshold:';

            % Create lcmv_thEditField
            app.lcmv_thEditField = uieditfield(app.activ_MethodsPanel, 'numeric');
            app.lcmv_thEditField.Limits = [0 1.7321];
            app.lcmv_thEditField.FontSize = 14;
            app.lcmv_thEditField.Position = [528 155 53 22];
            app.lcmv_thEditField.Value = 1;

            % Create gamma1EditField_2Label
            app.gamma1EditField_2Label = uilabel(app.activ_MethodsPanel);
            app.gamma1EditField_2Label.HorizontalAlignment = 'right';
            app.gamma1EditField_2Label.FontSize = 14;
            app.gamma1EditField_2Label.Position = [455 115 63 22];
            app.gamma1EditField_2Label.Text = 'gamma1:';

            % Create lcmv_gamma1EditField
            app.lcmv_gamma1EditField = uieditfield(app.activ_MethodsPanel, 'numeric');
            app.lcmv_gamma1EditField.Limits = [0 2];
            app.lcmv_gamma1EditField.FontSize = 14;
            app.lcmv_gamma1EditField.Position = [529 115 53 22];

            % Create gamma2EditField_2Label
            app.gamma2EditField_2Label = uilabel(app.activ_MethodsPanel);
            app.gamma2EditField_2Label.HorizontalAlignment = 'right';
            app.gamma2EditField_2Label.FontSize = 14;
            app.gamma2EditField_2Label.Position = [455 76 63 22];
            app.gamma2EditField_2Label.Text = 'gamma2:';

            % Create lcmv_gamma2EditField
            app.lcmv_gamma2EditField = uieditfield(app.activ_MethodsPanel, 'numeric');
            app.lcmv_gamma2EditField.Limits = [0 2];
            app.lcmv_gamma2EditField.FontSize = 14;
            app.lcmv_gamma2EditField.Position = [529 76 53 22];
            app.lcmv_gamma2EditField.Value = 2;

            % Create deltagammaEditField_2Label
            app.deltagammaEditField_2Label = uilabel(app.activ_MethodsPanel);
            app.deltagammaEditField_2Label.HorizontalAlignment = 'right';
            app.deltagammaEditField_2Label.FontSize = 14;
            app.deltagammaEditField_2Label.Position = [428 37 90 22];
            app.deltagammaEditField_2Label.Text = 'delta gamma:';

            % Create lcmv_deltagammaEditField
            app.lcmv_deltagammaEditField = uieditfield(app.activ_MethodsPanel, 'numeric');
            app.lcmv_deltagammaEditField.Limits = [0.1 1];
            app.lcmv_deltagammaEditField.FontSize = 14;
            app.lcmv_deltagammaEditField.Position = [529 37 53 22];
            app.lcmv_deltagammaEditField.Value = 0.1;

            % Create activ_OtherparamsPanel
            app.activ_OtherparamsPanel = uipanel(app.ActivationTab);
            app.activ_OtherparamsPanel.Title = 'Other params';
            app.activ_OtherparamsPanel.FontWeight = 'bold';
            app.activ_OtherparamsPanel.FontSize = 14;
            app.activ_OtherparamsPanel.Position = [10 51 675 175];

            % Create IsNeighLabel
            app.IsNeighLabel = uilabel(app.activ_OtherparamsPanel);
            app.IsNeighLabel.FontSize = 14;
            app.IsNeighLabel.Tooltip = {'default <<false>>. no neighbor structure <<false>>, Laplacian neighbor structure  <<true>>.'};
            app.IsNeighLabel.Position = [420 117 56 22];
            app.IsNeighLabel.Text = 'IsNeigh:';

            % Create IsParcelLabel
            app.IsParcelLabel = uilabel(app.activ_OtherparamsPanel);
            app.IsParcelLabel.FontSize = 14;
            app.IsParcelLabel.Tooltip = {'default <<true>>. Structured sparsity with smoothness of the responses within areas <<false>> (no smoothness) <<true>> (parcel smoothness)'};
            app.IsParcelLabel.Position = [239 118 59 22];
            app.IsParcelLabel.Text = 'IsParcel:';

            % Create IsCurvLabel
            app.IsCurvLabel = uilabel(app.activ_OtherparamsPanel);
            app.IsCurvLabel.FontSize = 14;
            app.IsCurvLabel.Tooltip = {'Depth compensation factor for the Lead Field of giri and sulci <<false>> (no compensation) <<true>> (compensation)'};
            app.IsCurvLabel.Position = [58 118 49 22];
            app.IsCurvLabel.Text = 'IsCurv:';

            % Create activ_IsNeighCheckBox
            app.activ_IsNeighCheckBox = uicheckbox(app.activ_OtherparamsPanel);
            app.activ_IsNeighCheckBox.Tooltip = {'default <<false>>. no neighbor structure <<false>>, Laplacian neighbor structure  <<true>>.'};
            app.activ_IsNeighCheckBox.Text = '';
            app.activ_IsNeighCheckBox.FontSize = 14;
            app.activ_IsNeighCheckBox.Position = [484 115 26 22];

            % Create activ_IsParcelCheckBox
            app.activ_IsParcelCheckBox = uicheckbox(app.activ_OtherparamsPanel);
            app.activ_IsParcelCheckBox.Tooltip = {'default <<true>>. Structured sparsity with smoothness of the responses within areas <<false>> (no smoothness) <<true>> (parcel smoothness)'};
            app.activ_IsParcelCheckBox.Text = '';
            app.activ_IsParcelCheckBox.FontSize = 14;
            app.activ_IsParcelCheckBox.Position = [305 118 26 22];

            % Create activ_IsCurvCheckBox
            app.activ_IsCurvCheckBox = uicheckbox(app.activ_OtherparamsPanel);
            app.activ_IsCurvCheckBox.Tooltip = {'Depth compensation factor for the Lead Field of giri and sulci <<false>> (no compensation) <<true>> (compensation)'};
            app.activ_IsCurvCheckBox.Text = '';
            app.activ_IsCurvCheckBox.FontSize = 14;
            app.activ_IsCurvCheckBox.Position = [113 116 39 24];

            % Create activ_aGiriEditField
            app.activ_aGiriEditField = uieditfield(app.activ_OtherparamsPanel, 'numeric');
            app.activ_aGiriEditField.Limits = [1 5];
            app.activ_aGiriEditField.FontSize = 14;
            app.activ_aGiriEditField.Tooltip = {'default <<5>>. Baseline of giri curvature factor'};
            app.activ_aGiriEditField.Position = [110 74 55 22];
            app.activ_aGiriEditField.Value = 5;

            % Create aGiriLabel
            app.aGiriLabel = uilabel(app.activ_OtherparamsPanel);
            app.aGiriLabel.HorizontalAlignment = 'right';
            app.aGiriLabel.FontSize = 14;
            app.aGiriLabel.Position = [57 74 38 22];
            app.aGiriLabel.Text = 'aGiri:';

            % Create activ_aSulcEditField
            app.activ_aSulcEditField = uieditfield(app.activ_OtherparamsPanel, 'numeric');
            app.activ_aSulcEditField.Limits = [1 5];
            app.activ_aSulcEditField.FontSize = 14;
            app.activ_aSulcEditField.Tooltip = {'default <<5>>. Baseline of sulci curvature factor. usually is set to be higher than aGiri'};
            app.activ_aSulcEditField.Position = [110 18 55 22];
            app.activ_aSulcEditField.Value = 5;

            % Create aSulcEditFieldLabel
            app.aSulcEditFieldLabel = uilabel(app.activ_OtherparamsPanel);
            app.aSulcEditFieldLabel.HorizontalAlignment = 'right';
            app.aSulcEditFieldLabel.FontSize = 14;
            app.aSulcEditFieldLabel.Position = [53 18 44 22];
            app.aSulcEditFieldLabel.Text = 'aSulc:';

            % Create activ_bGiriEditField
            app.activ_bGiriEditField = uieditfield(app.activ_OtherparamsPanel, 'numeric');
            app.activ_bGiriEditField.Limits = [0 5];
            app.activ_bGiriEditField.FontSize = 14;
            app.activ_bGiriEditField.Tooltip = {'default <<3>>. Scale of giri curvature factor'};
            app.activ_bGiriEditField.Position = [304 74 60 22];
            app.activ_bGiriEditField.Value = 3;

            % Create bGiriEditFieldLabel
            app.bGiriEditFieldLabel = uilabel(app.activ_OtherparamsPanel);
            app.bGiriEditFieldLabel.HorizontalAlignment = 'right';
            app.bGiriEditFieldLabel.FontSize = 14;
            app.bGiriEditFieldLabel.Position = [253 74 38 22];
            app.bGiriEditFieldLabel.Text = 'bGiri:';

            % Create activ_bSulcEditField
            app.activ_bSulcEditField = uieditfield(app.activ_OtherparamsPanel, 'numeric');
            app.activ_bSulcEditField.Limits = [1 5];
            app.activ_bSulcEditField.FontSize = 14;
            app.activ_bSulcEditField.Tooltip = {'default <<3>>. Scale of sulci curvature factor'};
            app.activ_bSulcEditField.Position = [304 18 60 22];
            app.activ_bSulcEditField.Value = 3;

            % Create bSulcEditFieldLabel
            app.bSulcEditFieldLabel = uilabel(app.activ_OtherparamsPanel);
            app.bSulcEditFieldLabel.HorizontalAlignment = 'right';
            app.bSulcEditFieldLabel.FontSize = 14;
            app.bSulcEditFieldLabel.Position = [249 18 44 22];
            app.bSulcEditFieldLabel.Text = 'bSulc:';

            % Create activ_regLaplacianEditField
            app.activ_regLaplacianEditField = uieditfield(app.activ_OtherparamsPanel, 'numeric');
            app.activ_regLaplacianEditField.Limits = [0.1 1];
            app.activ_regLaplacianEditField.FontSize = 14;
            app.activ_regLaplacianEditField.Position = [484 17 60 22];
            app.activ_regLaplacianEditField.Value = 0.1;

            % Create regLaplacianEditFieldLabel
            app.regLaplacianEditFieldLabel = uilabel(app.activ_OtherparamsPanel);
            app.regLaplacianEditFieldLabel.HorizontalAlignment = 'right';
            app.regLaplacianEditFieldLabel.FontSize = 14;
            app.regLaplacianEditFieldLabel.Position = [382 17 89 22];
            app.regLaplacianEditFieldLabel.Text = 'regLaplacian:';

            % Create activ_fieldSlider
            app.activ_fieldSlider = uislider(app.activ_OtherparamsPanel);
            app.activ_fieldSlider.Limits = [1 3];
            app.activ_fieldSlider.MajorTicks = [1 2 3];
            app.activ_fieldSlider.MinorTicks = [1 2 3];
            app.activ_fieldSlider.Tooltip = {'default <<2>>. Use projected Lead Field or 3D profile with prior of rotational invariance <<1>> (projected Lead Field) <<2>> (2D field isotropy) <<3>> (3D field isotropy)'};
            app.activ_fieldSlider.FontSize = 14;
            app.activ_fieldSlider.Position = [483 88 92 3];
            app.activ_fieldSlider.Value = 2;

            % Create FieldLabel
            app.FieldLabel = uilabel(app.activ_OtherparamsPanel);
            app.FieldLabel.HorizontalAlignment = 'right';
            app.FieldLabel.FontSize = 14;
            app.FieldLabel.Position = [431 78 39 22];
            app.FieldLabel.Text = 'Field:';

            % Create ConnectivityTab
            app.ConnectivityTab = uitab(app.TabGroup);
            app.ConnectivityTab.Tooltip = {'Parameters in the third level of analysis or HIGGS connectivity. if the second level is targeted then the first level will be performed in source screening mode'};
            app.ConnectivityTab.Title = 'Connectivity';

            % Create conn_MethodsPanel
            app.conn_MethodsPanel = uipanel(app.ConnectivityTab);
            app.conn_MethodsPanel.Title = 'Methods';
            app.conn_MethodsPanel.FontWeight = 'bold';
            app.conn_MethodsPanel.FontSize = 14;
            app.conn_MethodsPanel.Position = [9 372 674 106];

            % Create HIGGSLabel
            app.HIGGSLabel = uilabel(app.conn_MethodsPanel);
            app.HIGGSLabel.FontSize = 14;
            app.HIGGSLabel.FontWeight = 'bold';
            app.HIGGSLabel.Position = [120 53 55 22];
            app.HIGGSLabel.Text = 'HIGGS:';

            % Create CheckBox_higgs
            app.CheckBox_higgs = uicheckbox(app.conn_MethodsPanel);
            app.CheckBox_higgs.Text = '';
            app.CheckBox_higgs.FontSize = 14;
            app.CheckBox_higgs.Position = [180 53 26 22];

            % Create ThresholdEditFieldLabel_2
            app.ThresholdEditFieldLabel_2 = uilabel(app.conn_MethodsPanel);
            app.ThresholdEditFieldLabel_2.HorizontalAlignment = 'right';
            app.ThresholdEditFieldLabel_2.FontSize = 14;
            app.ThresholdEditFieldLabel_2.Position = [98 17 71 22];
            app.ThresholdEditFieldLabel_2.Text = 'Threshold:';

            % Create higgs_thEditField
            app.higgs_thEditField = uieditfield(app.conn_MethodsPanel, 'numeric');
            app.higgs_thEditField.Limits = [0 1.7321];
            app.higgs_thEditField.FontSize = 14;
            app.higgs_thEditField.Position = [177 17 54 22];
            app.higgs_thEditField.Value = 1;

            % Create hg_LASSOLabel
            app.hg_LASSOLabel = uilabel(app.conn_MethodsPanel);
            app.hg_LASSOLabel.FontSize = 14;
            app.hg_LASSOLabel.FontWeight = 'bold';
            app.hg_LASSOLabel.Position = [322 53 83 22];
            app.hg_LASSOLabel.Text = 'hg_LASSO:';

            % Create CheckBox_hg_lasso
            app.CheckBox_hg_lasso = uicheckbox(app.conn_MethodsPanel);
            app.CheckBox_hg_lasso.Text = '';
            app.CheckBox_hg_lasso.FontSize = 14;
            app.CheckBox_hg_lasso.Position = [406 53 26 22];

            % Create ThresholdEditField_2Label_2
            app.ThresholdEditField_2Label_2 = uilabel(app.conn_MethodsPanel);
            app.ThresholdEditField_2Label_2.HorizontalAlignment = 'right';
            app.ThresholdEditField_2Label_2.FontSize = 14;
            app.ThresholdEditField_2Label_2.Position = [329 15 71 22];
            app.ThresholdEditField_2Label_2.Text = 'Threshold:';

            % Create hg_lasso_thEditField
            app.hg_lasso_thEditField = uieditfield(app.conn_MethodsPanel, 'numeric');
            app.hg_lasso_thEditField.Limits = [0 1.7321];
            app.hg_lasso_thEditField.FontSize = 14;
            app.hg_lasso_thEditField.Position = [405 15 53 22];
            app.hg_lasso_thEditField.Value = 1;

            % Create conn_OtherparamsPanel
            app.conn_OtherparamsPanel = uipanel(app.ConnectivityTab);
            app.conn_OtherparamsPanel.Title = 'Other params';
            app.conn_OtherparamsPanel.FontWeight = 'bold';
            app.conn_OtherparamsPanel.FontSize = 14;
            app.conn_OtherparamsPanel.Position = [10 49 673 310];

            % Create conn_IsCurvCheckBox
            app.conn_IsCurvCheckBox = uicheckbox(app.conn_OtherparamsPanel);
            app.conn_IsCurvCheckBox.Tooltip = {'default <<true>>. Depth compensation factor for the Lead Field of giri and sulci <<false>> (no compensation) <<true>> (compensation)'};
            app.conn_IsCurvCheckBox.Text = '';
            app.conn_IsCurvCheckBox.FontSize = 14;
            app.conn_IsCurvCheckBox.Position = [98 251 26 22];

            % Create conn_IsNeighCheckBox
            app.conn_IsNeighCheckBox = uicheckbox(app.conn_OtherparamsPanel);
            app.conn_IsNeighCheckBox.Tooltip = {'default <<false>>. no neighbor structure <<false>>. Laplacian neighbor structure  <<true>>.'};
            app.conn_IsNeighCheckBox.Text = '';
            app.conn_IsNeighCheckBox.FontSize = 14;
            app.conn_IsNeighCheckBox.Position = [314 252 26 22];

            % Create aGiriEditField_2Label
            app.aGiriEditField_2Label = uilabel(app.conn_OtherparamsPanel);
            app.aGiriEditField_2Label.HorizontalAlignment = 'right';
            app.aGiriEditField_2Label.FontSize = 14;
            app.aGiriEditField_2Label.Position = [40 207 38 22];
            app.aGiriEditField_2Label.Text = 'aGiri:';

            % Create conn_aGiriEditField
            app.conn_aGiriEditField = uieditfield(app.conn_OtherparamsPanel, 'numeric');
            app.conn_aGiriEditField.Limits = [1 5];
            app.conn_aGiriEditField.FontSize = 14;
            app.conn_aGiriEditField.Tooltip = {'default <<5>>. Baseline of giri curvature factor'};
            app.conn_aGiriEditField.Position = [97 207 68 22];
            app.conn_aGiriEditField.Value = 5;

            % Create aSulcEditField_2Label
            app.aSulcEditField_2Label = uilabel(app.conn_OtherparamsPanel);
            app.aSulcEditField_2Label.HorizontalAlignment = 'right';
            app.aSulcEditField_2Label.FontSize = 14;
            app.aSulcEditField_2Label.Position = [36 158 44 22];
            app.aSulcEditField_2Label.Text = 'aSulc:';

            % Create conn_aSulcEditField
            app.conn_aSulcEditField = uieditfield(app.conn_OtherparamsPanel, 'numeric');
            app.conn_aSulcEditField.Limits = [1 5];
            app.conn_aSulcEditField.FontSize = 14;
            app.conn_aSulcEditField.Tooltip = {'default <<5>>. Baseline of sulci curvature factor. usually is set to be higher than aGiri'};
            app.conn_aSulcEditField.Position = [95 158 70 22];
            app.conn_aSulcEditField.Value = 5;

            % Create bGiriEditField_2Label
            app.bGiriEditField_2Label = uilabel(app.conn_OtherparamsPanel);
            app.bGiriEditField_2Label.HorizontalAlignment = 'right';
            app.bGiriEditField_2Label.FontSize = 14;
            app.bGiriEditField_2Label.Position = [262 207 38 22];
            app.bGiriEditField_2Label.Text = 'bGiri:';

            % Create conn_bGiriEditField
            app.conn_bGiriEditField = uieditfield(app.conn_OtherparamsPanel, 'numeric');
            app.conn_bGiriEditField.Limits = [0 5];
            app.conn_bGiriEditField.FontSize = 14;
            app.conn_bGiriEditField.Tooltip = {'default <<3>>. Scale of giri curvature factor'};
            app.conn_bGiriEditField.Position = [315 207 68 22];
            app.conn_bGiriEditField.Value = 3;

            % Create bSulcEditField_2Label
            app.bSulcEditField_2Label = uilabel(app.conn_OtherparamsPanel);
            app.bSulcEditField_2Label.HorizontalAlignment = 'right';
            app.bSulcEditField_2Label.FontSize = 14;
            app.bSulcEditField_2Label.Position = [258 158 44 22];
            app.bSulcEditField_2Label.Text = 'bSulc:';

            % Create conn_bSulcEditField
            app.conn_bSulcEditField = uieditfield(app.conn_OtherparamsPanel, 'numeric');
            app.conn_bSulcEditField.Limits = [1 5];
            app.conn_bSulcEditField.FontSize = 14;
            app.conn_bSulcEditField.Tooltip = {'default <<3>>. Scale of sulci curvature factor'};
            app.conn_bSulcEditField.Position = [315 158 68 22];
            app.conn_bSulcEditField.Value = 3;

            % Create regLaplacianEditField_2Label
            app.regLaplacianEditField_2Label = uilabel(app.conn_OtherparamsPanel);
            app.regLaplacianEditField_2Label.HorizontalAlignment = 'right';
            app.regLaplacianEditField_2Label.FontSize = 14;
            app.regLaplacianEditField_2Label.Position = [214 23 89 22];
            app.regLaplacianEditField_2Label.Text = 'regLaplacian:';

            % Create conn_regLaplacianEditField
            app.conn_regLaplacianEditField = uieditfield(app.conn_OtherparamsPanel, 'numeric');
            app.conn_regLaplacianEditField.Limits = [0.1 1];
            app.conn_regLaplacianEditField.FontSize = 14;
            app.conn_regLaplacianEditField.Position = [316 23 67 22];
            app.conn_regLaplacianEditField.Value = 0.1;

            % Create MaxiterouterEditFieldLabel
            app.MaxiterouterEditFieldLabel = uilabel(app.conn_OtherparamsPanel);
            app.MaxiterouterEditFieldLabel.HorizontalAlignment = 'right';
            app.MaxiterouterEditFieldLabel.FontSize = 14;
            app.MaxiterouterEditFieldLabel.Position = [428 207 91 22];
            app.MaxiterouterEditFieldLabel.Text = 'Maxiter outer:';

            % Create conn_MaxiterouterEditField
            app.conn_MaxiterouterEditField = uieditfield(app.conn_OtherparamsPanel, 'numeric');
            app.conn_MaxiterouterEditField.Limits = [20 80];
            app.conn_MaxiterouterEditField.FontSize = 14;
            app.conn_MaxiterouterEditField.Tooltip = {'default <<60>>. Maximum number of iterations for the HIGGS EM outer loop'};
            app.conn_MaxiterouterEditField.Position = [534 207 71 22];
            app.conn_MaxiterouterEditField.Value = 60;

            % Create MaxiterinnerEditFieldLabel
            app.MaxiterinnerEditFieldLabel = uilabel(app.conn_OtherparamsPanel);
            app.MaxiterinnerEditFieldLabel.HorizontalAlignment = 'right';
            app.MaxiterinnerEditFieldLabel.FontSize = 14;
            app.MaxiterinnerEditFieldLabel.Position = [429 159 90 22];
            app.MaxiterinnerEditFieldLabel.Text = 'Maxiter inner:';

            % Create conn_MaxiterinnerEditField
            app.conn_MaxiterinnerEditField = uieditfield(app.conn_OtherparamsPanel, 'numeric');
            app.conn_MaxiterinnerEditField.Limits = [20 40];
            app.conn_MaxiterinnerEditField.FontSize = 14;
            app.conn_MaxiterinnerEditField.Tooltip = {'default <<30>>. Maximum number of iterations for the HIGGS hGLASSO inner loop'};
            app.conn_MaxiterinnerEditField.Position = [534 159 71 22];
            app.conn_MaxiterinnerEditField.Value = 30;

            % Create ntryLabel
            app.ntryLabel = uilabel(app.conn_OtherparamsPanel);
            app.ntryLabel.HorizontalAlignment = 'right';
            app.ntryLabel.FontSize = 14;
            app.ntryLabel.Position = [488 52 32 22];
            app.ntryLabel.Text = 'ntry:';

            % Create conn_ntrySlider
            app.conn_ntrySlider = uislider(app.conn_OtherparamsPanel);
            app.conn_ntrySlider.Limits = [0 5];
            app.conn_ntrySlider.MajorTicks = [0 1 2 3 4 5];
            app.conn_ntrySlider.MinorTicks = [0 1 2 3 4 5];
            app.conn_ntrySlider.Tooltip = {'default <<3>>. Number of iterations forward in eveidence HIGGS eveidence optimality prediction to evaluate the performance of the hGLASSO Rayleigh threshold at every iteration'};
            app.conn_ntrySlider.FontSize = 14;
            app.conn_ntrySlider.Position = [534 69 116 3];

            % Create axiLabel
            app.axiLabel = uilabel(app.conn_OtherparamsPanel);
            app.axiLabel.HorizontalAlignment = 'right';
            app.axiLabel.FontSize = 14;
            app.axiLabel.Position = [53 110 27 22];
            app.axiLabel.Text = 'axi:';

            % Create conn_axiEditField
            app.conn_axiEditField = uieditfield(app.conn_OtherparamsPanel, 'numeric');
            app.conn_axiEditField.Limits = [0.1 1];
            app.conn_axiEditField.FontSize = 14;
            app.conn_axiEditField.Tooltip = {'default <<1E-1>>. Instrumental noise inferior limit that can be defined according to experimental information, <<1E-1>> represents 10%'};
            app.conn_axiEditField.Position = [94 110 71 22];
            app.conn_axiEditField.Value = 0.1;

            % Create conn_prewarmingCheckBox
            app.conn_prewarmingCheckBox = uicheckbox(app.conn_OtherparamsPanel);
            app.conn_prewarmingCheckBox.Tooltip = {'default <<1>>. Perform prewarming of the HIGGS EM by means of a two-step methdod (first activation and then connectivity). <<0>> (no prewarming) <<1>> (does prewarming)'};
            app.conn_prewarmingCheckBox.Text = '';
            app.conn_prewarmingCheckBox.FontSize = 14;
            app.conn_prewarmingCheckBox.Position = [534 251 26 22];

            % Create penaltySliderLabel
            app.penaltySliderLabel = uilabel(app.conn_OtherparamsPanel);
            app.penaltySliderLabel.HorizontalAlignment = 'right';
            app.penaltySliderLabel.FontSize = 14;
            app.penaltySliderLabel.Position = [466 112 54 22];
            app.penaltySliderLabel.Text = 'penalty:';

            % Create conn_penaltySlider
            app.conn_penaltySlider = uislider(app.conn_OtherparamsPanel);
            app.conn_penaltySlider.Limits = [0 2];
            app.conn_penaltySlider.MajorTicks = [0 1 2];
            app.conn_penaltySlider.MinorTicks = [0 1 2];
            app.conn_penaltySlider.Tooltip = {'default <<2>>. Penalization model for HIGGS <<0>> (naive or vareta) <<1>> (hermitian graphical LASSO or hGLASSO) <<2>> (hermitian graphical Ridge or hgRidge). if set to <<1>> the computational cost is higher but with full statistical guaranties. <<2>> does not ofer full guaranties but can be very similar to <<1>>'};
            app.conn_penaltySlider.FontSize = 14;
            app.conn_penaltySlider.Position = [533 128 81 3];

            % Create rth1Label
            app.rth1Label = uilabel(app.conn_OtherparamsPanel);
            app.rth1Label.HorizontalAlignment = 'right';
            app.rth1Label.FontSize = 14;
            app.rth1Label.Position = [47 66 33 22];
            app.rth1Label.Text = 'rth1:';

            % Create conn_rth1EditField
            app.conn_rth1EditField = uieditfield(app.conn_OtherparamsPanel, 'numeric');
            app.conn_rth1EditField.Limits = [0.1 1];
            app.conn_rth1EditField.FontSize = 14;
            app.conn_rth1EditField.Tooltip = {'default <<0.7>>. Inferior Rayleigh threshold <<0.7>> is the maximum of the Rayleigh distribution'};
            app.conn_rth1EditField.Position = [93 66 72 22];
            app.conn_rth1EditField.Value = 0.7;

            % Create rth2Label
            app.rth2Label = uilabel(app.conn_OtherparamsPanel);
            app.rth2Label.HorizontalAlignment = 'right';
            app.rth2Label.FontSize = 14;
            app.rth2Label.Position = [269 61 33 22];
            app.rth2Label.Text = 'rth2:';

            % Create conn_rth2EditField
            app.conn_rth2EditField = uieditfield(app.conn_OtherparamsPanel, 'numeric');
            app.conn_rth2EditField.Limits = [0.1 5];
            app.conn_rth2EditField.FontSize = 14;
            app.conn_rth2EditField.Tooltip = {'default <<3.16>>. Superior Rayleigh threshold <<3.16>> is the 99% percentile of the Rayleigh distribution'};
            app.conn_rth2EditField.Position = [315 61 68 22];
            app.conn_rth2EditField.Value = 3.16;

            % Create eigregLabel
            app.eigregLabel = uilabel(app.conn_OtherparamsPanel);
            app.eigregLabel.HorizontalAlignment = 'right';
            app.eigregLabel.FontSize = 14;
            app.eigregLabel.Position = [32 25 48 22];
            app.eigregLabel.Text = 'eigreg:';

            % Create conn_eigregEditField
            app.conn_eigregEditField = uieditfield(app.conn_OtherparamsPanel, 'numeric');
            app.conn_eigregEditField.Limits = [0.0001 0.1];
            app.conn_eigregEditField.FontSize = 14;
            app.conn_eigregEditField.Position = [93 25 72 22];
            app.conn_eigregEditField.Value = 0.0001;

            % Create FieldLabel_2
            app.FieldLabel_2 = uilabel(app.conn_OtherparamsPanel);
            app.FieldLabel_2.HorizontalAlignment = 'right';
            app.FieldLabel_2.FontSize = 14;
            app.FieldLabel_2.Position = [263 114 39 22];
            app.FieldLabel_2.Text = 'Field:';

            % Create conn_fieldSlider
            app.conn_fieldSlider = uislider(app.conn_OtherparamsPanel);
            app.conn_fieldSlider.Limits = [1 3];
            app.conn_fieldSlider.MajorTicks = [1 2 3];
            app.conn_fieldSlider.MinorTicks = [1 2 3];
            app.conn_fieldSlider.Tooltip = {'default <<2>>. Use projected Lead Field or 3D profile with prior of rotational invariance <<1>> (projected Lead Field) <<2>> (2D field isotropy)'};
            app.conn_fieldSlider.FontSize = 14;
            app.conn_fieldSlider.Position = [315 131 84 3];
            app.conn_fieldSlider.Value = 2;

            % Create IsCurvLabel_2
            app.IsCurvLabel_2 = uilabel(app.conn_OtherparamsPanel);
            app.IsCurvLabel_2.FontSize = 14;
            app.IsCurvLabel_2.Tooltip = {'default <<true>>. Depth compensation factor for the Lead Field of giri and sulci <<false>> (no compensation) <<true>> (compensation)'};
            app.IsCurvLabel_2.Position = [35 253 49 22];
            app.IsCurvLabel_2.Text = 'IsCurv:';

            % Create IsNeighLabel_2
            app.IsNeighLabel_2 = uilabel(app.conn_OtherparamsPanel);
            app.IsNeighLabel_2.FontSize = 14;
            app.IsNeighLabel_2.Tooltip = {'default <<false>>. no neighbor structure <<false>>. Laplacian neighbor structure  <<true>>.'};
            app.IsNeighLabel_2.Position = [251 253 56 22];
            app.IsNeighLabel_2.Text = 'IsNeigh:';

            % Create prewarmingLabel
            app.prewarmingLabel = uilabel(app.conn_OtherparamsPanel);
            app.prewarmingLabel.FontSize = 14;
            app.prewarmingLabel.Tooltip = {'default <<1>>. Perform prewarming of the HIGGS EM by means of a two-step methdod (first activation and then connectivity). <<0>> (no prewarming) <<1>> (does prewarming)'};
            app.prewarmingLabel.Position = [443 253 82 22];
            app.prewarmingLabel.Text = 'prewarming:';

            % Create RunButton
            app.RunButton = uibutton(app.BCVARETAToolboxv10UIFigure, 'push');
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @RunButtonPushed, true);
            app.RunButton.Icon = fullfile(pathToMLAPP, 'images', 'run.jpg');
            app.RunButton.HorizontalAlignment = 'left';
            app.RunButton.FontSize = 14;
            app.RunButton.Position = [606 15 78 25];
            app.RunButton.Text = 'Run';

            % Create CancelButton
            app.CancelButton = uibutton(app.BCVARETAToolboxv10UIFigure, 'push');
            app.CancelButton.ButtonPushedFcn = createCallbackFcn(app, @CancelButtonPushed, true);
            app.CancelButton.Icon = fullfile(pathToMLAPP, 'images', 'cancel.jpg');
            app.CancelButton.HorizontalAlignment = 'left';
            app.CancelButton.FontSize = 14;
            app.CancelButton.Position = [507 15 90 25];
            app.CancelButton.Text = 'Cancel';

            % Show the figure after all components are created
            app.BCVARETAToolboxv10UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = BC_VARETA

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.BCVARETAToolboxv10UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.BCVARETAToolboxv10UIFigure)
        end
    end
end