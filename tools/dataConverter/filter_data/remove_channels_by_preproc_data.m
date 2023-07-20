function [channel_layout,leadfield] = remove_channels_by_preproc_data(labels,channel_layout,leadfield)

from = 1;
limit = length(channel_layout.Channel);

if(nargin<4)
    while(from <= limit)
        pos = find(strcmpi(channel_layout.Channel(from).Name, labels), 1);
        if (isempty(pos))
            channel_layout.Channel(from)=[];
            leadfield(from,:)=[];
            limit = limit - 1;
        else
            from = from + 1;
        end
    end
else
    while(from <= limit)
        if(isnan(leadfield(from,1)))
            channel_layout.Channel(from)=[];
            leadfield(from,:)=[];
            limit = limit - 1;
        else
            from = from + 1;
        end
    end
end
end

