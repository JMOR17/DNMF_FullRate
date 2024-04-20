function [cROIs, Cs, coherence, skew, sz, patchID, tElapsed, dimensions] = mcb_DNMF_v5(path_to_video, options)
    % [cROIs, Cs, coherence, skew, sz, tElapsed] = mcb_DNMF(path_to_video, options)
    
    % Load in Video
    fileName = path_to_video;
    Y = bigread2(fileName,1);
    if(size(Y,3)==1)
        Y = imread_big(fileName);
    end
    % Set saturated pixels to 0
    Y(Y>options.maxVal) = 0;
    
    tStart = tic;
    [szA, szB, ~] = size(Y);
    
    if(nargin<2)
        options = defaultOptions_mcbDNMF();
    end
    
    minSkew = options.minSkew;
    sizeRange = options.sizeRange;
        
    
    % Detect ROIs in patches
    [cROIs0, Cs0, coherence0, skew0, sz0, patchID0] = DNMF_General5(Y, options);
    
    % Light screen for skewness and ROI size
    valid = skew0>minSkew & isbetween(sz0,sizeRange(1),sizeRange(2));
    cROIs = cROIs0(:,valid);
    Cs = Cs0(valid,:);
    coherence = coherence0(valid);
    skew = skew0(valid);
    sz = sz0(valid);
    patchID = patchID0(valid);
    dimensions = [szA szB];
    
    %% No merging. Leave that to the user later
    tElapsed = toc(tStart);
    fprintf('Done!\n');
end