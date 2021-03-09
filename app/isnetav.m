function x=isnetav
disp('-->> Checking internet status');
[~,b]=dos('ping -n 1 www.github.io');
n=strfind(b,'Lost');
n1=b(n+7);
if(n1=='0')
    x=1;
    disp('-->> Internet connection Ok');
else
    x=0;
    disp('-->> We could not get you online');
    disp('-->> There migth be either a problem with the internet connection or matlab proxy configuration');
end
end
