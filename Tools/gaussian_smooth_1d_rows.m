function [sm_data] = gaussian_smooth_1d_rows(data,sig,windur)

    sm_data = NaN(size(data));

    if nargin < 3        
        windur = 3; %default
    end

    for i=1:size(data,1)
        sm_data(i,:) = gaussian_smooth_1d(data(i,:),sig,windur);
        
    end

end