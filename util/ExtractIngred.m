function [recipes_binary_subsamp,ingreds_subsamp,sorted_ingred_freq]=...
  ExtractIngred(Ningred_used,ingred_freq,recipes_binary,ingreds)
% Extracts the Ningred_used most common ingredients from a recipe data set
%
% Written by: DJ Strouse
% Last updated: Aug 19, 2013 by DJ Strouse
%
% INPUTS
% Ningred_used [=] positive integer = number of ingredients to extract
% ingred_freq [=] Ningred-length vector of positive integers = number of
%     recipes that each ingredient appears in
% recipes_binary [=] Nrec X Ningred binary matrix = indicates appearance of
%     ingredients in recipes
% ingreds [=] Ningred-element array = names of ingredients
%
% OUTPUTS
% recipes_binary_subsamp
% ingreds_subsamp
% sorted_ingred_freq

[sorted_ingred_freq,i] = sort(ingred_freq,'descend');
to_use = i(1:Ningred_used); clear i;
recipes_binary_subsamp = recipes_binary(:,to_use)';
ingreds_subsamp = ingreds(to_use); clear to_use;
disp(sprintf('Using the following %i ingredients:',Ningred_used));
disp(ingreds_subsamp');

end

