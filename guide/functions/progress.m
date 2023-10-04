function progress
    fig = uifigure;
    d = uiprogressdlg(fig,'Title','Please Wait',...
        'Message','Opening the application');
    pause(.5)

    % Perform calculations
    % ...
    d.Value = .33; 
    d.Message = 'Loading your data';
    pause(1)

    % Perform calculations
    % ...
    d.Value = .67;
    d.Message = 'Processing the data';
    pause(1)

    % Finish calculations
    % ...
    d.Value = 1;
    d.Message = 'Finishing';
    pause(1)

    % Close dialog box
    close(d)
end
