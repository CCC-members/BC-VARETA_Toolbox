function processImage(doc, reportData, holeId, varargin)
%This function adds an image to our report
import mlreportgen.dom.*

for i=1:length(varargin)
    eval([inputname(i+3) '= varargin{i};']);
end

%We access the reportData structure dynamically
fileName = reportData.(holeId);
%A new Image objected is created and the image is sized 10x10cm
img = Image(fileName);

clasif = split(holeId,'_');
if(isequal(clasif{2},'bandtopography'))
    holeId = 'bandtopography';
end
if(isequal(clasif{2},'source'))
    holeId = 'source';
end

switch holeId
    case 'Image_activation'
        img.Width  = '6.33in';
        img.Height = '3.95in';
    case 'Image_headmodel_front'
        img.Width  = '3in';
        img.Height = '2.95in';
    case 'Image_headmodel_top'
        img.Width  = '3in';
        img.Height = '2.95in';
    case 'Image_headmodel_right'
        img.Width  = '3in';
        img.Height = '2.95in';
    case 'Image_headmodel_left'
        img.Width  = '3in';
        img.Height = '2.95in';
        case 'Image_sourcemodel_top'
        img.Width  = '3in';
        img.Height = '2.95in';
    case 'Image_sourcemodel_bottom'
        img.Width  = '3in';
        img.Height = '2.95in';
    case 'Image_sourcemodel_right'
        img.Width  = '3in';
        img.Height = '2.95in';
    case 'Image_sourcemodel_left'
        img.Width  = '3in';
        img.Height = '2.95in';
        case 'Image_leadfield_top'
        img.Width  = '3in';
        img.Height = '2.95in';
    case 'Image_leadfield_front'
        img.Width  = '3in';
        img.Height = '2.95in';
    case 'Image_leadfield_right'
        img.Width  = '3in';
        img.Height = '2.95in';
    case 'Image_leadfield_left'
        img.Width  = '3in';
        img.Height = '2.95in';
    case 'Image_channels_front'
        img.Width  = '3in';
        img.Height = '2.95in';
    case 'Image_channels_back'
        img.Width  = '3in';
        img.Height = '2.95in';
    case 'Image_channels_right'
        img.Width  = '3in';
        img.Height = '2.95in';
    case 'Image_channels_left'
        img.Width  = '3in';
        img.Height = '2.95in';
    case 'Image_sensorPSD'
        img.Width  = '6.47in';
        img.Height = '4.0in';
    case 'bandtopography'
        if(nfreqs==4)
            img.Width  = '3.2in';
            img.Height = '3.2in';
        elseif(nfreqs==5 || nfreqs==6)
            img.Width  = '2.2in';
            img.Height = '2.2in';
        else
            img.Width  = '2in';
            img.Height = '2in';
        end   
    case 'source'
        if(nfreqs==4)
            img.Width  = '3.2in';
            img.Height = '3.2in';
        elseif(nfreqs==5)
            img.Width  = '2in';
            img.Height = '1.6in';
        else
            img.Width  = '2in';
            img.Height = '2in';
        end   
    case 'Image_alphathetaration1'
        img.Width  = '3in';
        img.Height = '3in'; 
    case 'Image_alphathetaration2'
        img.Width  = '3in';
        img.Height = '3in';
    case 'Image_alphapeakpositionPSD'
        img.Width  = '3.45in';
        img.Height = '2.8in';
    case 'Image_alphapeakposition'
        img.Width  = '2.6in';
        img.Height = '3in';    
    case 'Image_normativedist'
        img.Width  = '6.3in';
        img.Height = '3.95in';
    case 'Image_deltacoronal'
        img.Width  = '2.7in';
        img.Height = '2.7in';
    case 'Image_deltalateral'
        img.Width  = '3.26in';
        img.Height = '2.8in';
    case 'Image_thetacoronal'
        img.Width  = '2.7in';
        img.Height = '2.7in';
    case 'Image_thetalateral'
        img.Width  = '3.26in';
        img.Height = '2.8in';
    case 'Image_alphacoronal'
        img.Width  = '2.7in';
        img.Height = '2.7in';
    case 'Image_alphalateral'
        img.Width  = '3.26in';
        img.Height = '2.8in';
    case 'Image_betacoronal'
        img.Width  = '2.7in';
        img.Height = '2.7in';
    case 'Image_betalateral'
        img.Width  = '3.26in';
        img.Height = '2.8in';
    case 'Image_activityconnectivity_part_01'
        img.Width  = '6.25in';
        img.Height = '3.95in';
    case 'Image_activityconnectivity_part_02'
        img.Width  = '6.25in';
        img.Height = '3.95in';
    case 'Image_activityconnectivity_part_03'
        img.Width  = '6.25in';
        img.Height = '3.95in'; 
    otherwise
        img.Width  = '8in';
        img.Height = '8in';
end
%We create a new paragraph, based on the AR_Image style and append the
%image to this paragraph. This is done because this style has a center
%alignment, and so the image is centered too.
p = Paragraph( '', 'AR_Image' );
append(p, img);
append(doc, p);

%After the image we add another paragraph of AR_Caption style. This
%style counts the number of images and adds "Figure X: " on front of our
%text
% p = Paragraph('This is the MATLAB logo.', 'AR_Caption');
% append(doc, p);
end


