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
    disp('Platform not supported');
end
if(n1=='5')
    x=1;
    disp('-->> Internet connection Ok');
else
    x=0;
    cprintf('red','-->> We could not get you online'); fprintf('\n');
    cprintf('red','-->> There migth be either a problem with the internet connection or matlab proxy configuration'); fprintf('\n');
    cprintf('red','-->> We could not get the BC-VARETA latest version right now.');fprintf('\n');
    cprintf('red','-->> Please check your internet configuration and try to update again.');fprintf('\n');
    disp('=================================================================');
end
end
