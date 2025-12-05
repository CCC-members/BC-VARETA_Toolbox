function EEG = remove_eeg_channels_by_labels(user_labels, EEG)
data        = EEG.data;
labels      = {EEG.chanlocs.labels};
from        = 1;
limit       = size(data,1);
clean_labels = {size(user_labels,1),1};
for i=1:length(user_labels)
    clean_labels{i} = strrep(user_labels{i},' ','');
end

while(from <= limit)
    pos = find(strcmpi(labels{from}, clean_labels), 1);
    if (isempty(pos))
        data(from,:)    = [];
        labels(from)    = [];
        limit           = limit - 1;        
    else
        from = from + 1;
    end
end
EEG.data                    = data;
rej_indms                   = length(labels)+1:length(EEG.chanlocs);
EEG.chanlocs(rej_indms)     = [];
[EEG.chanlocs.labels]       = labels{:};
EEG.nbchan                  = length(labels);
end

