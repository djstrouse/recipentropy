function [ output ] = rectify( input )
% rectify sets all negative elements of a scalar, vector, or matrix input
% to zero and returns the resulting (rectified) object

pos_entries = input>0;
output = pos_entries .* input;

end

