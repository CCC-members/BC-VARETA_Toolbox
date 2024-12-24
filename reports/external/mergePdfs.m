% Merges the pdf-Documents in the input cell array fileNames into one
% single pdf-Document with file name outputFile

function mergePdfs(fileNames, outputFile)

memSet = org.apache.pdfbox.io.MemoryUsageSetting.setupMainMemoryOnly();
merger = org.apache.pdfbox.multipdf.PDFMergerUtility;

cellfun(@(f) merger.addSource(f), fileNames)

merger.setDestinationFileName(outputFile)
merger.mergeDocuments(memSet)
