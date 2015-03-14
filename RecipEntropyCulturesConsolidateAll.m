function [] = RecipEntropyCulturesConsolidateAll(culture,Ningred_used,consolidate_orders,Nsubsamps)
% Consolidates maxent data for Ningred_used most common ingredients
% in recipes from a specified culture
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
% Ningred_used [=] positive integer = number of spices to use
% consolidate_orders [=] boolean = indicates whether or not to loop over
%     runs of RecipEntropyConsolidate.m first
% Nsubsamps [=] vector of non-negative integers = Nsubsamp for each order
%
% OUTPUTS
% none

% init
if ~exist('consolidate_orders','var')
  consolidate_orders = true; % defaults to true
end
if ~exist('Nsubsamps','var')
  Nsubsamps = [25 25 25];
end
if length(Nsubsamps)==1
  Nsubsamps = Nsubsamps*ones(3,1);
end
if length(Nsubsamps)~=3
  error('Nsubsamps should have 3 elements, one for each order!')
end
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
if consolidate_orders
  orders = 0:2;
  RecipEntropyCulturesConsolidate(culture,Ningred_used,orders,Nsubsamps);
end

% load recipe data set to extract its data
load('data/recipes.mat');
Ningred = size(recipes_binary,2);
Nrec = size(recipes_binary,1);
[~,recipes_binary,Nrec_used]=...
  ExtractCulturalRecipes(culture,recipes,recipes_binary);
subsamp_size = round(.8*Nrec_used);
[recipes_binary_subsamp,ingreds_subsamp] =...
  ExtractIngred(Ningred_used,ingred_freq,recipes_binary,ingreds);

% load freqs and maxents
load(['data/cultures/',culture_str,sprintf('/%iingred/order%i/large.mat',Ningred_used,0)]);
load(['data/cultures/',culture_str,sprintf('/%iingred/order%i/small.mat',Ningred_used,0)]);
load(['data/cultures/',culture_str,sprintf('/%iingred/order%i/large.mat',Ningred_used,1)]);
load(['data/cultures/',culture_str,sprintf('/%iingred/order%i/small.mat',Ningred_used,1)]);
maxent1_small = maxent_small; clear maxent_small;
maxent1_large = maxent_large; clear maxent_large;
load(['data/cultures/',culture_str,sprintf('/%iingred/order%i/large.mat',Ningred_used,2)]);
load(['data/cultures/',culture_str,sprintf('/%iingred/order%i/small.mat',Ningred_used,2)]);
maxent2_small = maxent_small; clear maxent_small;
maxent2_large = maxent_large; clear maxent_large;

% check for consistency in Ningred, Nrec, Nsubsamp, and subsamp_size


% compare maxents
maxent1v2 = CompareMaxEnt12(...
  maxent1_large,maxent2_large,freq_large,...
  'first-order maxent','second-order maxent','raw freqs');

% save
save(...
  ['data/cultures/',culture_str,sprintf('/%iingred/large.mat',Ningred_used)],...
  'maxent1_large',...
  'maxent2_large',...
  'freq_large',...
  'maxent1v2',...
  'culture',...
  'culture_str',...
  'Ningred',...
  'Nrec',...
  'Nrec_used',...
  'Nsubsamps',...
  'subsamp_size',...
  'Ningred_used',...
  'ingreds_subsamp',...
  'recipes_binary_subsamp',...
  '-v7.3');
save(...
  ['data/cultures/',culture_str,sprintf('/%iingred/small.mat',Ningred_used)],...
  'maxent1_small',...
  'maxent2_small',...
  'freq_small',...
  'maxent1v2',...
  'culture',...
  'culture_str',...
  'Ningred',...
  'Nrec',...
  'Nrec_used',...
  'Nsubsamps',...
  'subsamp_size',...
  'Ningred_used',...
  '-v7.3');

end