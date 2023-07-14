function [properties] = define_frequency_bands(properties)

%--------------------frequency bands-----------------------------------
if( properties.general_params.run_frequency_bin.value)
    freq_resol = properties.sensor_params.freq_resol.value;
    frequency_name                      = {};
    frequency_bin                       = {};
    frequency_run                       = {};
    band_start                          = {};
    band_end                            = {};
    str_band                            = {};
    if(getGlobalGuimode())
        process_bin_waitbar =  waitbar(0,strcat('Computing the frequency''s bin...'));
    end
    frequency_bands = properties.sensor_params.frequencies;
    pos = 1;
    for h=1:size(frequency_bands,1)
        if(getGlobalGuimode())
            waitbar(h/size(frequency_bands,1),process_bin_waitbar,strcat('Computing the frequency''s bin...'));
        end
        band = frequency_bands(h);
        if(band.run)
            pointer                     = band.f_start;
            frequency_name{pos}         = band.name;
            frequency_bin{pos}          = pointer;
            frequency_run{pos}          = band.run;
            band_start{pos}             = band.f_start;
            band_end{pos}               = band.f_end;
            str_band{pos}               = strcat(band.name,"_",string(pointer),"Hz");
            pos                         = pos + 1;

            while band.f_end > round( pointer + freq_resol , 4 )
                pointer                 = pointer + freq_resol;
                frequency_name{pos}     = band.name;
                frequency_bin{pos}      = pointer;
                frequency_run{pos}      = band.run;
                band_start{pos}         = band.f_start;
                band_end{pos}           = band.f_end;
                str_band{pos}           = strcat(band.name,"_",string(pointer),"Hz");
                pos                     = pos + 1;
            end
            if(pointer < band.f_end)
                pointer                 = pointer + freq_resol;
                frequency_name{pos}     = band.name;
                frequency_bin{pos}      = band.f_end;
                frequency_run{pos}      = band.run;
                band_start{pos}         = band.f_start;
                band_end{pos}           = band.f_end;
                str_band{pos}           = strcat(band.name,'_',string(pointer),'Hz');
                pos                     = pos + 1;
            end
        end
    end
    if(getGlobalGuimode())
        delete(process_bin_waitbar);
    end
    properties.sensor_params.frequencies = struct('name',frequency_name,'f_bin',frequency_bin,'f_start',band_start,'f_end',band_end,'run',frequency_run,'str_band',str_band);
else
    frequency_name                      = {};
    frequency_run                       = {};
    band_start                          = {};
    band_end                            = {};
    str_band                            = {};
    frequency_bands                     = properties.sensor_params.frequencies;
    pos                                 = 1;
    for h=1:size(frequency_bands,1)
        band                            = frequency_bands(h);
        if(band.run)
            frequency_name{pos}         = band.name;
            frequency_run{pos}          = band.run;
            band_start{pos}             = band.f_start;
            band_end{pos}               = band.f_end;
            str_band{pos}               = strcat( band.name,'_',string(band.f_start),'Hz_',string(band.f_end),'Hz');
            pos                         = pos + 1;
        end
    end
    properties.sensor_params.frequencies = struct('name',frequency_name,'f_start',band_start,'f_end',band_end,'run',frequency_run,'str_band',str_band);
end
end


