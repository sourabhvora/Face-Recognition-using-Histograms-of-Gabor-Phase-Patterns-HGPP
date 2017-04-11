function score = histogram_intersection( h, s, num_bins )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    h = reshape(h,num_bins,[]);
    s = reshape(s,num_bins,[]);
    
    score = sum(sum( min(h,s) ));
    %score = score/sum(h);

end

