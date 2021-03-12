function FLAG=test_gpu(count,memorySize)
if nargin<2
    memorySize=4e9;
end

if gpuDeviceCount>0
    if nargin==1 && ~isempty(count)
        gpuInfo=gpuDevice(count);
    else
        try
            gpuInfo=gpuDevice();
        catch
            fprintf('gpuDevice() has error, skipped.\n');
            gpuInfo=gpuDevice(1);
        end
    end
    if gpuInfo.AvailableMemory>memorySize
        FLAG=true;
    else
        FLAG=false;
    end
else
    FLAG=false;
end
%     FLAG=false;

end
