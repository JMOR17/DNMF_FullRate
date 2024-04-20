function [Cs_sm, temp] = prep_FR(Cs, SM_AMOUNT, DETREND_FRAMES)
    
    if(nargin<3)
        DETREND_FRAMES = 900;
    end
    if(nargin<2)
        SM_AMOUNT = 3;
    end
    
    DETREND_DIMENSION = 2;
    Cs_sm = gaussian_smooth_1d_rows(Cs,SM_AMOUNT);
    [~, temp] = j_detrend2b(Cs_sm, DETREND_FRAMES, DETREND_DIMENSION, false, true);
    Cs_sm = Cs_sm-temp;


end