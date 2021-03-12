classdef freqresol_maxfreq_samplfreq_guide < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        SpectralpropertiesUIFigure     matlab.ui.Figure
        Panel                          matlab.ui.container.Panel
        SamplingfrequencySpinnerLabel  matlab.ui.control.Label
        SamplingfrequencySpinner       matlab.ui.control.Spinner
        FrequencyresolutionSpinnerLabel  matlab.ui.control.Label
        FrequencyresolutionSpinner     matlab.ui.control.Spinner
        MaximumfrequencySpinnerLabel   matlab.ui.control.Label
        MaximumfrequencySpinner        matlab.ui.control.Spinner
        WindowingorderSpinnerLabel     matlab.ui.control.Label
        WindowingorderSpinner          matlab.ui.control.Spinner
        CancelButton                   matlab.ui.control.Button
        OkButton                       matlab.ui.control.Button
    end

    
    properties (Access = private)
        Property % Description
    end
    
    properties (Access = public)
        canceled
        frequency_resolution % Resolution frequency
        sampling_frequency % Samplin frequency
        max_frequency % Max frequency
        
    end
    

    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            
        end

        % Button pushed function: CancelButton
        function CancelButtonPushed(app, event)
            app.canceled = true;
            uiresume(app.SpectralpropertiesUIFigure);
            close();
        end

        % Button pushed function: OkButton
        function OkButtonPushed(app, event)
            app.canceled = false;
            bcv_properties = jsondecode(fileread(fullfile('bcv_properties','bcv_properties.json')));
                       
            bcv_properties.freq_resol.value = app.FrequencyresolutionSpinner.Value;
            bcv_properties.samp_freq.value = app.SamplingfrequencySpinner.Value;
            bcv_properties.max_freq.value = app.MaximumfrequencySpinner.Value;
            bcv_properties.win_order.value = app.WindowingorderSpinner.Value;
                                
            saveJSON(bcv_properties,fullfile('bcv_properties','bcv_properties.json'));
            
            uiresume(app.SpectralpropertiesUIFigure);
            
        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create SpectralpropertiesUIFigure
            app.SpectralpropertiesUIFigure = uifigure;
            app.SpectralpropertiesUIFigure.Position = [100 100 318 284];
            app.SpectralpropertiesUIFigure.Name = 'Spectral properties';

            % Create Panel
            app.Panel = uipanel(app.SpectralpropertiesUIFigure);
            app.Panel.Position = [32 69 258 194];

            % Create SamplingfrequencySpinnerLabel
            app.SamplingfrequencySpinnerLabel = uilabel(app.Panel);
            app.SamplingfrequencySpinnerLabel.HorizontalAlignment = 'right';
            app.SamplingfrequencySpinnerLabel.Position = [21 113 115 22];
            app.SamplingfrequencySpinnerLabel.Text = 'Sampling frequency:';

            % Create SamplingfrequencySpinner
            app.SamplingfrequencySpinner = uispinner(app.Panel);
            app.SamplingfrequencySpinner.Limits = [0 1000];
            app.SamplingfrequencySpinner.Position = [151 113 71 22];
            app.SamplingfrequencySpinner.Value = 200;

            % Create FrequencyresolutionSpinnerLabel
            app.FrequencyresolutionSpinnerLabel = uilabel(app.Panel);
            app.FrequencyresolutionSpinnerLabel.HorizontalAlignment = 'right';
            app.FrequencyresolutionSpinnerLabel.Position = [15 155 121 22];
            app.FrequencyresolutionSpinnerLabel.Text = 'Frequency resolution:';

            % Create FrequencyresolutionSpinner
            app.FrequencyresolutionSpinner = uispinner(app.Panel);
            app.FrequencyresolutionSpinner.Step = 0.1;
            app.FrequencyresolutionSpinner.Position = [151 155 71 22];
            app.FrequencyresolutionSpinner.Value = 0.5;

            % Create MaximumfrequencySpinnerLabel
            app.MaximumfrequencySpinnerLabel = uilabel(app.Panel);
            app.MaximumfrequencySpinnerLabel.HorizontalAlignment = 'right';
            app.MaximumfrequencySpinnerLabel.Position = [19 66 117 22];
            app.MaximumfrequencySpinnerLabel.Text = 'Maximum frequency:';

            % Create MaximumfrequencySpinner
            app.MaximumfrequencySpinner = uispinner(app.Panel);
            app.MaximumfrequencySpinner.Step = 0.1;
            app.MaximumfrequencySpinner.Limits = [0 Inf];
            app.MaximumfrequencySpinner.Position = [151 66 71 22];
            app.MaximumfrequencySpinner.Value = 32;

            % Create WindowingorderSpinnerLabel
            app.WindowingorderSpinnerLabel = uilabel(app.Panel);
            app.WindowingorderSpinnerLabel.HorizontalAlignment = 'right';
            app.WindowingorderSpinnerLabel.Position = [38 24 99 22];
            app.WindowingorderSpinnerLabel.Text = 'Windowing order:';

            % Create WindowingorderSpinner
            app.WindowingorderSpinner = uispinner(app.Panel);
            app.WindowingorderSpinner.Limits = [1 10];
            app.WindowingorderSpinner.Position = [152 24 71 22];
            app.WindowingorderSpinner.Value = 3;

            % Create CancelButton
            app.CancelButton = uibutton(app.SpectralpropertiesUIFigure, 'push');
            app.CancelButton.ButtonPushedFcn = createCallbackFcn(app, @CancelButtonPushed, true);
            app.CancelButton.Position = [211 32 79 22];
            app.CancelButton.Text = 'Cancel';

            % Create OkButton
            app.OkButton = uibutton(app.SpectralpropertiesUIFigure, 'push');
            app.OkButton.ButtonPushedFcn = createCallbackFcn(app, @OkButtonPushed, true);
            app.OkButton.Position = [116 32 85 22];
            app.OkButton.Text = 'Ok';
        end
    end

    methods (Access = public)

        % Construct app
        function app = freqresol_maxfreq_samplfreq_guide

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.SpectralpropertiesUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.SpectralpropertiesUIFigure)
        end
    end
end