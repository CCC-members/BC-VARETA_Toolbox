classdef frequency_bands_guide < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                   matlab.ui.Figure
        SelectfrequencysbandPanel  matlab.ui.container.Panel
        CheckBox_Delta             matlab.ui.control.CheckBox
        CheckBox_theta             matlab.ui.control.CheckBox
        CheckBox_alpha             matlab.ui.control.CheckBox
        CheckBox_beta              matlab.ui.control.CheckBox
        FromLabel                  matlab.ui.control.Label
        HztoLabel                  matlab.ui.control.Label
        HzdeltaLabel               matlab.ui.control.Label
        Spinner_delta_start        matlab.ui.control.Spinner
        Spinner_delta_end          matlab.ui.control.Spinner
        FromLabel_2                matlab.ui.control.Label
        HzthetaLabel               matlab.ui.control.Label
        FromLabel_3                matlab.ui.control.Label
        HzalphaLabel               matlab.ui.control.Label
        FromLabel_4                matlab.ui.control.Label
        HzbetaLabel                matlab.ui.control.Label
        HztoLabel_2                matlab.ui.control.Label
        Spinner_theta_start        matlab.ui.control.Spinner
        Spinner_theta_end          matlab.ui.control.Spinner
        HztoLabel_3                matlab.ui.control.Label
        Spinner_alpha_start        matlab.ui.control.Spinner
        Spinner_alpha_end          matlab.ui.control.Spinner
        HztoLabel_4                matlab.ui.control.Label
        Spinner_beta_start         matlab.ui.control.Spinner
        Spinner_beta_end           matlab.ui.control.Spinner
        AgreeButton                matlab.ui.control.Button
        CancelButton               matlab.ui.control.Button
        RunmodeButtonGroup         matlab.ui.container.ButtonGroup
        ByfrequencybandButton      matlab.ui.control.RadioButton
        FrequencybinButton         matlab.ui.control.RadioButton
    end

    
    properties (Access = public)
        frequencies % Description
        Property2 % Description
        canceled
        frequency_bin
    end
    
    methods (Access = private)
        
        function results = func(app)
            
        end
    end
    

    methods (Access = private)

        % Button pushed function: AgreeButton
        function AgreeButtonPushed(app, event)
            app.canceled = false;
            
            bcv_properties = jsondecode(fileread(fullfile('bcv_properties','bcv_properties.json')));
            
            bcv_properties.run_frequency_bin.value = app.FrequencybinButton.Value;
            
            for i = 1 : length(bcv_properties.frequencies)
                switch bcv_properties.frequencies(i).name
                    case 'delta'
                        if(app.CheckBox_Delta.Value)
                            if(app.Spinner_delta_start.Value >= app.Spinner_delta_end.Value)
                                uiwait(msgbox('The start value to the delta band must be smaller than end value','Error','modal'));
                                return;
                            end
                            bcv_properties.frequencies(i).f_start = app.Spinner_delta_start.Value;
                            bcv_properties.frequencies(i).f_end = app.Spinner_delta_end.Value;
                        end
                        bcv_properties.frequencies(i).run = app.CheckBox_Delta.Value;
                    case 'theta'
                        if(app.CheckBox_theta.Value)
                            if(app.Spinner_theta_start.Value >= app.Spinner_theta_end.Value)
                                uiwait(msgbox('The start value to the theta band must be smaller than end value','Error','modal'));
                                return;
                            end
                            bcv_properties.frequencies(i).f_start = app.Spinner_theta_start.Value;
                            bcv_properties.frequencies(i).f_end = app.Spinner_theta_end.Value;
                        end
                        bcv_properties.frequencies(i).run = app.CheckBox_theta.Value;
                    case 'alpha'
                        if(app.CheckBox_alpha.Value)
                            if(app.Spinner_alpha_start.Value >= app.Spinner_alpha_end.Value)
                                uiwait(msgbox('The start value to the alpha band must be smaller than end value','Error','modal'));
                                return;
                            end
                            bcv_properties.frequencies(i).f_start = app.Spinner_alpha_start.Value;
                            bcv_properties.frequencies(i).f_end = app.Spinner_alpha_end.Value;
                        end
                        bcv_properties.frequencies(i).run = app.CheckBox_alpha.Value;
                    case 'beta'
                        if(app.CheckBox_beta.Value)
                            if(app.Spinner_beta_start.Value >= app.Spinner_beta_end.Value)
                                uiwait(msgbox('The start value to the beta band must be smaller than end value','Error','modal'));
                                return;
                            end
                            bcv_properties.frequencies(i).f_start = app.Spinner_beta_start.Value;
                            bcv_properties.frequencies(i).f_end = app.Spinner_beta_end.Value;
                        end
                        bcv_properties.frequencies(i).run = app.CheckBox_beta.Value;
                end
            end
            
            saveJSON(bcv_properties,fullfile('bcv_properties','bcv_properties.json'));
            
            app.canceled = false;
            uiresume(app.UIFigure);
            
        end

        % Button pushed function: CancelButton
        function CancelButtonPushed(app, event)
            app.frequencies = [];
            app.canceled = true;
            uiresume(app.UIFigure);
            close();
            
            
            
        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure
            app.UIFigure = uifigure;
            app.UIFigure.Position = [100 100 368 356];
            app.UIFigure.Name = 'UI Figure';

            % Create SelectfrequencysbandPanel
            app.SelectfrequencysbandPanel = uipanel(app.UIFigure);
            app.SelectfrequencysbandPanel.Title = 'Select frequency''s band';
            app.SelectfrequencysbandPanel.Position = [29 62 309 191];

            % Create CheckBox_Delta
            app.CheckBox_Delta = uicheckbox(app.SelectfrequencysbandPanel);
            app.CheckBox_Delta.Text = '';
            app.CheckBox_Delta.Position = [26 134 25 22];
            app.CheckBox_Delta.Value = true;

            % Create CheckBox_theta
            app.CheckBox_theta = uicheckbox(app.SelectfrequencysbandPanel);
            app.CheckBox_theta.Text = '';
            app.CheckBox_theta.Position = [26 101 25 22];
            app.CheckBox_theta.Value = true;

            % Create CheckBox_alpha
            app.CheckBox_alpha = uicheckbox(app.SelectfrequencysbandPanel);
            app.CheckBox_alpha.Text = '';
            app.CheckBox_alpha.Position = [26 64 25 22];
            app.CheckBox_alpha.Value = true;

            % Create CheckBox_beta
            app.CheckBox_beta = uicheckbox(app.SelectfrequencysbandPanel);
            app.CheckBox_beta.Text = '';
            app.CheckBox_beta.Position = [26 30 25 22];
            app.CheckBox_beta.Value = true;

            % Create FromLabel
            app.FromLabel = uilabel(app.SelectfrequencysbandPanel);
            app.FromLabel.Position = [50 134 33 22];
            app.FromLabel.Text = 'From';

            % Create HztoLabel
            app.HztoLabel = uilabel(app.SelectfrequencysbandPanel);
            app.HztoLabel.Position = [140 134 40 22];
            app.HztoLabel.Text = ' Hz to ';

            % Create HzdeltaLabel
            app.HzdeltaLabel = uilabel(app.SelectfrequencysbandPanel);
            app.HzdeltaLabel.Position = [237 134 58 22];
            app.HzdeltaLabel.Text = 'Hz (delta)';

            % Create Spinner_delta_start
            app.Spinner_delta_start = uispinner(app.SelectfrequencysbandPanel);
            app.Spinner_delta_start.Step = 0.1;
            app.Spinner_delta_start.Limits = [0 50];
            app.Spinner_delta_start.Position = [84 134 57 22];
            app.Spinner_delta_start.Value = 0.1;

            % Create Spinner_delta_end
            app.Spinner_delta_end = uispinner(app.SelectfrequencysbandPanel);
            app.Spinner_delta_end.Step = 0.1;
            app.Spinner_delta_end.Limits = [0.1 50];
            app.Spinner_delta_end.Position = [174 134 53 22];
            app.Spinner_delta_end.Value = 4;

            % Create FromLabel_2
            app.FromLabel_2 = uilabel(app.SelectfrequencysbandPanel);
            app.FromLabel_2.Position = [51 101 33 22];
            app.FromLabel_2.Text = 'From';

            % Create HzthetaLabel
            app.HzthetaLabel = uilabel(app.SelectfrequencysbandPanel);
            app.HzthetaLabel.Position = [236 101 58 22];
            app.HzthetaLabel.Text = 'Hz (theta)';

            % Create FromLabel_3
            app.FromLabel_3 = uilabel(app.SelectfrequencysbandPanel);
            app.FromLabel_3.Position = [52 64 33 22];
            app.FromLabel_3.Text = 'From';

            % Create HzalphaLabel
            app.HzalphaLabel = uilabel(app.SelectfrequencysbandPanel);
            app.HzalphaLabel.Position = [235 64 61 22];
            app.HzalphaLabel.Text = 'Hz (alpha)';

            % Create FromLabel_4
            app.FromLabel_4 = uilabel(app.SelectfrequencysbandPanel);
            app.FromLabel_4.Position = [52 30 33 22];
            app.FromLabel_4.Text = 'From';

            % Create HzbetaLabel
            app.HzbetaLabel = uilabel(app.SelectfrequencysbandPanel);
            app.HzbetaLabel.Position = [235 30 55 22];
            app.HzbetaLabel.Text = 'Hz (beta)';

            % Create HztoLabel_2
            app.HztoLabel_2 = uilabel(app.SelectfrequencysbandPanel);
            app.HztoLabel_2.Position = [140 100 40 22];
            app.HztoLabel_2.Text = ' Hz to ';

            % Create Spinner_theta_start
            app.Spinner_theta_start = uispinner(app.SelectfrequencysbandPanel);
            app.Spinner_theta_start.Step = 0.1;
            app.Spinner_theta_start.Limits = [0 50];
            app.Spinner_theta_start.Position = [84 100 57 22];
            app.Spinner_theta_start.Value = 4;

            % Create Spinner_theta_end
            app.Spinner_theta_end = uispinner(app.SelectfrequencysbandPanel);
            app.Spinner_theta_end.Step = 0.1;
            app.Spinner_theta_end.Limits = [0.1 50];
            app.Spinner_theta_end.Position = [174 100 53 22];
            app.Spinner_theta_end.Value = 7;

            % Create HztoLabel_3
            app.HztoLabel_3 = uilabel(app.SelectfrequencysbandPanel);
            app.HztoLabel_3.Position = [140 64 40 22];
            app.HztoLabel_3.Text = ' Hz to ';

            % Create Spinner_alpha_start
            app.Spinner_alpha_start = uispinner(app.SelectfrequencysbandPanel);
            app.Spinner_alpha_start.Step = 0.1;
            app.Spinner_alpha_start.Limits = [0 50];
            app.Spinner_alpha_start.Position = [84 64 57 22];
            app.Spinner_alpha_start.Value = 7;

            % Create Spinner_alpha_end
            app.Spinner_alpha_end = uispinner(app.SelectfrequencysbandPanel);
            app.Spinner_alpha_end.Step = 0.1;
            app.Spinner_alpha_end.Limits = [0.1 50];
            app.Spinner_alpha_end.Position = [174 64 53 22];
            app.Spinner_alpha_end.Value = 14;

            % Create HztoLabel_4
            app.HztoLabel_4 = uilabel(app.SelectfrequencysbandPanel);
            app.HztoLabel_4.Position = [140 30 40 22];
            app.HztoLabel_4.Text = ' Hz to ';

            % Create Spinner_beta_start
            app.Spinner_beta_start = uispinner(app.SelectfrequencysbandPanel);
            app.Spinner_beta_start.Step = 0.1;
            app.Spinner_beta_start.Limits = [0 50];
            app.Spinner_beta_start.Position = [84 30 57 22];
            app.Spinner_beta_start.Value = 14;

            % Create Spinner_beta_end
            app.Spinner_beta_end = uispinner(app.SelectfrequencysbandPanel);
            app.Spinner_beta_end.Step = 0.1;
            app.Spinner_beta_end.Limits = [0.1 50];
            app.Spinner_beta_end.Position = [174 30 53 22];
            app.Spinner_beta_end.Value = 31;

            % Create AgreeButton
            app.AgreeButton = uibutton(app.UIFigure, 'push');
            app.AgreeButton.ButtonPushedFcn = createCallbackFcn(app, @AgreeButtonPushed, true);
            app.AgreeButton.Position = [118 19 100 22];
            app.AgreeButton.Text = 'Agree';

            % Create CancelButton
            app.CancelButton = uibutton(app.UIFigure, 'push');
            app.CancelButton.ButtonPushedFcn = createCallbackFcn(app, @CancelButtonPushed, true);
            app.CancelButton.Position = [238 19 100 22];
            app.CancelButton.Text = 'Cancel';

            % Create RunmodeButtonGroup
            app.RunmodeButtonGroup = uibuttongroup(app.UIFigure);
            app.RunmodeButtonGroup.Title = 'Run mode';
            app.RunmodeButtonGroup.Position = [29 273 309 60];

            % Create ByfrequencybandButton
            app.ByfrequencybandButton = uiradiobutton(app.RunmodeButtonGroup);
            app.ByfrequencybandButton.Text = 'By frequency band ';
            app.ByfrequencybandButton.Position = [11 14 125 22];
            app.ByfrequencybandButton.Value = true;

            % Create FrequencybinButton
            app.FrequencybinButton = uiradiobutton(app.RunmodeButtonGroup);
            app.FrequencybinButton.Text = 'Frequency bin';
            app.FrequencybinButton.Position = [188 14 98 22];
        end
    end

    methods (Access = public)

        % Construct app
        function app = frequency_bands_guide

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end