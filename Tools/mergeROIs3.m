function [cROIs, Cs] = mergeROIs3(Cs0, cROIs0, temporalCorrThr, spatialOverlapThr, fullSize, normalize_option, activity_option, SM_AMOUNT, DETREND_FRAMES)
   

%     fprintf('Merging ROIs...');

    Cs1 = prep_FR(Cs0, SM_AMOUNT, DETREND_FRAMES);
    % Z-score to compute temporal correlation
    
    zCs = zscore(Cs1,[],2);
    
    % Compute spatial overlap
    cROIs_binary = double(cROIs0>0);
    spatialOverlap_0 = cROIs_binary'*cROIs_binary;
    sizes = sum(cROIs_binary);
    [sz1, sz2] = meshgrid(sizes);
    temp = zeros(size(sz1,1),size(sz1,2),2);
    temp(:,:,1) = sz1;
    temp(:,:,2) = sz2;
    minSize = min(temp,[],3);

    spatialOverlap = spatialOverlap_0./minSize;
    if(activity_option)
        temporalCorr = zeros(size(spatialOverlap));
        for i1 = 1:size(Cs1b,1)-1
            for i2 = i1+1:size(Cs1b,1)
                valid = sum(zCs>2,1)>0;
                temp = corrcoef(zCs([i1 i2],valid)');
                temporalCorr(i1,i2) = temp(1,2);
                temporalCorr(i2,i1) = temp(1,2);
            end
        end
    else
        temporalCorr = zCs*zCs'/(size(zCs,2)-1);
    end
    temporalCorr(eye(size(Cs1,1))==1) = 1;
    % Find connected components of the graph
    linked = (spatialOverlap>spatialOverlapThr) & (temporalCorr>temporalCorrThr);    
    [~,connComp] = graphconncomp(sparse(linked));

    % Merge ROIs that are part of connected components
    cROIs = sparse(size(cROIs0,1), max(connComp));
    Cs = zeros(max(connComp), size(Cs1,2));
    coherence = NaN(max(connComp),1);
    skew = NaN(max(connComp),1);
    sz = NaN(max(connComp),1);
    for i_component=1:max(connComp)
        these = connComp==i_component; 
        if(sum(these)==1)                   % Connected components with a single member
            cROIs(:,i_component) = cROIs0(:,these);
            Cs(i_component,:) = Cs0(these,:);
        else                                % Connected components with multiple members
            rois = cROIs0(:,these);
            traces = Cs0(these,:);
            indices = sum(rois,2)>0;
            rois2 = rois(indices,:);
            yyy = rois2*traces;
            if(normalize_option)
                num = sum(rois2>0,2);
                yyy = yyy./num;
            end
            aaa = max(yyy,[],2);
            ccc = max(aaa\yyy,0);
            aaa2 = max(yyy*pinv(ccc),0);
            ccc2 = max(aaa2\yyy,0);
            n = sqrt(sum(aaa2.^2));
            aaa2 = aaa2/n;
            ccc2 = ccc2*n;
            cROIs(indices,i_component) = aaa2;
                        
            
            Cs(i_component,:) = ccc2;
        end
%         [ch, sk, z] = evaluateROIs(cROIs(:,i_component), Cs(i_component,:),fullSize);
%         coherence(i_component) = ch;
%         skew(i_component) = sk;
%         sz(i_component) = z;
    end
    
    
end