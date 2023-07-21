function [structuralPath,functionalPath,BC_V_Workspace] = download_testdata(testDataPath,BC_V_Workspace)

if(~isfolder(testDataPath) && ~isequal(testDataPath,fullfile(pwd,'data')))
    disp("The testData path is not a folder");
    return;
end
if(~isfolder(fullfile(pwd,'data')) && isequal(testDataPath,fullfile(pwd,'data')))
    mkdir(fullfile(pwd,'data'));
end
dafaultWorkspace = fullfile(char(java.lang.System.getProperty('user.home')),"BC-VARETA_Structure");
if(~isfolder(BC_V_Workspace) && ~isequal(BC_V_Workspace,dafaultWorkspace))
    disp("The Workspace path is not a folder");
    return;
end
if(~isfolder(BC_V_Workspace) )
    mkdir(BC_V_Workspace);
end

filename = fullfile(testDataPath,'TestData.zip');
if(~isfile(fullfile(testDataPath,'TestData.zip')))
    testDataURL = 'https://lstneuro-my.sharepoint.com/:u:/g/personal/ariosky_neuroinformatics-collaboratory_org/ETKFp83pAphAha99ka5xq3gBxH3b8CnC7KtWocKedF2RUA?e=srlmzT&download=1';    
    disp(strcat("-->> Downloading Test Data......."));
    matlab.net.http.HTTPOptions.VerifyServerName = false;
    options = weboptions('Timeout',Inf,'RequestMethod','auto');
    downladed_file = websave(filename,testDataURL,options);
end
%% Unzip lasted version
pause(1);
disp(strcat("-->> Unpacking test data..."));
unzip(filename,testDataPath);
pause(1);
%delete(filename);
disp("Test data is ready");

structuralPath  = fullfile(testDataPath,"Structural");
functionalPath  = fullfile(testDataPath,"Fuctional");




end

