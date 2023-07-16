function connected =connection_status()
connected = false;
try    
  options = weboptions('ContentType','json','Timeout',Inf,'RequestMethod','auto');  
  status = webread('https://www.github.com',options);        
  connected = true;
catch
    
end
end

