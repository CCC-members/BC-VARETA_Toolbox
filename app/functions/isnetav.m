function x=isnetav
disp('-->> Checking internet status');
[~,b]=system('ping -c 5 -W 15 www.github.com');
n=strfind(b,'received');
if ismac
    % Code to run on Mac platform
elseif isunix
   n1=b(n-2);
elseif ispc
   n1=b(n+7);
else
    disp('Platform not supported')
end
if(n1=='5')
    x=1;
    disp('-->> Internet connection Ok');
else
    x=0;
    disp('-->> We could not get you online');
    disp('-->> There migth be either a problem with the internet connection or matlab proxy configuration');
end
end
