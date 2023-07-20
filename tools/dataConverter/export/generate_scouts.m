function new_Scouts = generate_scouts(Cortex)
id_struc    = find(strcmp({Cortex.Atlas.Name} , "Structures") == 1);

Vertices        = num2cell(1:length(Cortex.Vertices));
Seed            = num2cell(1:length(Cortex.Vertices));
Color          = num2cell(rand(length(Cortex.Vertices),3),2)';
        
Label           = cellstr(string([1:length(Cortex.Vertices)]));
% Label           = char(Label);
Function        = string(ones(1,length(Cortex.Vertices)));
Function(:,1)   = "Mean";
Function        = cellstr(Function);
Region          = string(ones(1,length(Cortex.Vertices)));
Region(1,[1:fix(length(Cortex.Vertices)/2)]) = "LU";
Region(1,[fix(length(Cortex.Vertices)/2) + 1:length(Cortex.Vertices)]) = "RU";
Region           = cellstr(Region);
Handles          = string(zeros(1,length(Cortex.Vertices)));
new_Scouts      =  struct('Vertices',Vertices,'Seed',Seed,'Color',Color,'Label',Label,'Function',Function,'Region',Region);

% disp(strcat("   >> Surface: ",Cortex.Comment ));
% fprintf(1,'   >> Genering sCouts: %3d%%\n',0);
% for k=1:length(Cortex.Vertices)
%     sCout            = struct;
%     sCout.Vertices   = k;
%     sCout.Seed       = k;
%     sCout.Color      = [rand rand rand];
%     sCout.Label      = num2str(k);
%     sCout.Function   = 'Mean';
%     if(k<=length(Cortex.Vertices)/2)
%         sCout.Region = 'LU';
%     else
%         sCout.Region = 'RU';
%     end
%     sCout.Handles    = [];
%     new_Scouts(k)    = sCout;
%     fprintf(1,'\b\b\b\b%3.0f%%',(k)/(length(Cortex.Vertices))*100);
% end
end

