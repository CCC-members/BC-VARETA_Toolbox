function [structuralPath,functionalPath,outputPath] = download_testdata()
filename = strcat('data/TestData.zip');
if(~isfile(fullfile(pwd,"data",'TestData.zip')))
    testDataURL = 'https://lstneuro-my.sharepoint.com/:u:/g/personal/ariosky_neuroinformatics-collaboratory_org/ETKFp83pAphAha99ka5xq3gBxH3b8CnC7KtWocKedF2RUA?e=srlmzT&download=1';
    filename = strcat('data/TestData.zip');
    disp(strcat("-->> Downloading Test Data......."));
    matlab.net.http.HTTPOptions.VerifyServerName = false;
    options = weboptions('Timeout',Inf,'RequestMethod','auto');
    downladed_file = websave(filename,testDataURL,options);
end
%% Unzip lasted version
pause(1);
disp(strcat("-->> Unpacking test data..."));
exampleFiles = unzip(filename,'data/');
pause(1);
%delete(filename);
disp("Test data is ready");

structuralPath  = "data/Structural";
functionalPath  = "data/Fuctional";
outputPath      = fullfile(char(java.lang.System.getProperty('user.home')),"BC-VARETA_Structure");
if(~isfolder(outputPath))
mkdir(outputPath);
end

end

