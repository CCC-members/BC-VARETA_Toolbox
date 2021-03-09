function [folderpath,foldername] = unpackage_file(filename,unpackage_folder)
%UNPACKAGE_FILE Summary of this function goes here
%   Detailed explanation goes here
 f = dialog('Position',[300 300 250 80]);
 
 iconsClassName = 'com.mathworks.widgets.BusyAffordance$AffordanceSize';
 iconsSizeEnums = javaMethod('values',iconsClassName);
 SIZE_32x32 = iconsSizeEnums(2);  % (1) = 16x16,  (2) = 32x32
 jObj = com.mathworks.widgets.BusyAffordance(SIZE_32x32, 'Unpackaging file...');  % icon, label
 
 jObj.setPaintsWhenStopped(true);  % default = false
 jObj.useWhiteDots(false);         % default = false (true is good for dark backgrounds)
 javacomponent(jObj.getComponent, [50,10,150,80], f);
 jObj.start;
 pause(1);
 try
     folderpath = unzip(filename,unpackage_folder);
     [folderpath,foldername,ext] = fileparts(folderpath) ;
 catch
     delete(f);
     errordlg('Unpackage error!!!','Error');
     return;
 end
 jObj.stop;
 jObj.setBusyText('All done!');
 pause(2);
 delete(f);
 msgbox('Completed download!!!','Info');
                        
end

