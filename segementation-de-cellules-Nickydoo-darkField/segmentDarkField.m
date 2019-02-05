function [imOut, cellStats]=segmentDarkField(imInMt, imInAd, minObjectSize)
    mask = imbinarize(imInMt);
    imOut=bwareaopen(mask, minObjectSize);
    cellStats=struct2table(regionprops(imOut, 'Area', 'Centroid'));    
end
    