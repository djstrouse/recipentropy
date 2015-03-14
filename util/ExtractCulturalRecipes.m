function [recipes_cultural,recipes_binary_cultural,Nrec_used]=...
  ExtractCulturalRecipes(culture,recipes,recipes_binary)
% Extracts the recipes corresponding to a particular culture
%
% Written by: DJ Strouse
% Last updated: Aug 22, 2013 by DJ Strouse
%
% INPUTS
% culture [=] positive integer = ID of culture (ranked by number of recs)
%     1  = NorthAmerican
%     2  = SouthernEuropean
%     3  = LatinAmerican
%     4  = WesternEuropean
%     5  = EastAsian
%     6  = MiddleEastern
%     7  = SouthAsian
%     8  = SoutheastAsian
%     9  = EasternEuropean
%     10 = African
%     11 = NorthernEuropean
% recipes [=] struct = contains data about recipes
% recipes_binary [=] Nrec X Ningred binary matrix = indicates presence of
%     ingredients in recipes
%
% OUTPUTS
% recipes_cultural
% recipes_binary_cultural
% Nrec_used

% init
if culture==1
  culture_str = 'NorthAmerican';
elseif culture==2
  culture_str = 'SouthernEuropean';
elseif culture==3
  culture_str = 'LatinAmerican';
elseif culture==4
  culture_str = 'WesternEuropean';
elseif culture==5
  culture_str = 'EastAsian';
elseif culture==6
  culture_str = 'MiddleEastern';
elseif culture==7
  culture_str = 'SouthAsian';
elseif culture==8
  culture_str = 'SoutheastAsian';
elseif culture==9
  culture_str = 'EasternEuropean';
elseif culture==10
  culture_str = 'African';
elseif culture==11
  culture_str = 'NorthernEuropean';
else
  error('culture should be an integer between 1 and 11!')
end
Nrec = size(recipes_binary,1);

% extract recipes
culturemask = zeros(Nrec,1); % indicator function of culture identity
for r = 1:Nrec
  culturemask(r) = strcmp(recipes(r).region{1},culture_str);
end
clear r;
culturemask = culturemask==1; % convert to logical
recipes_cultural = recipes(culturemask);
recipes_binary_cultural = recipes_binary(culturemask,:);
Nrec_used = size(recipes_binary_cultural,1);

end

