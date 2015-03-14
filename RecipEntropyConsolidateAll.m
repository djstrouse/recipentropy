function [] = RecipEntropyConsolidateAll(Ningred_used,consolidate_orders,Nsubsamps)
% Consolidates maxent data for Ningred_used most common ingredients
%
% Written by: DJ Strouse
% Last updated: Aug 19, 2013 by DJ Strouse
%
% INPUTS
% Ningred_used [=] positive integer = number of ingredients to use
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
if consolidate_orders
  orders = 0:2;
  RecipEntropyConsolidate(Ningred_used,orders,Nsubsamps);
end

% load recipe data set to extract its data
load('data/recipes.mat');
Ningred = size(recipes_binary,2);
Nrec = size(recipes_binary,1);
subsamp_size = round(.8*Nrec);
[recipes_binary_subsamp,ingreds_subsamp] =...
  ExtractIngred(Ningred_used,ingred_freq,recipes_binary,ingreds);

% load freqs and maxents
load(sprintf('data/standard/%iingred/order%i/large.mat',Ningred_used,0));
load(sprintf('data/standard/%iingred/order%i/small.mat',Ningred_used,0));
load(sprintf('data/standard/%iingred/order%i/large.mat',Ningred_used,1));
load(sprintf('data/standard/%iingred/order%i/small.mat',Ningred_used,1));
maxent1_small = maxent_small; clear maxent_small;
maxent1_large = maxent_large; clear maxent_large;
load(sprintf('data/standard/%iingred/order%i/large.mat',Ningred_used,2));
load(sprintf('data/standard/%iingred/order%i/small.mat',Ningred_used,2));
maxent2_small = maxent_small; clear maxent_small;
maxent2_large = maxent_large; clear maxent_large;

% check for consistency in Ningred, Nrec, Nsubsamp, and subsamp_size


% compare maxents
maxent1v2 = CompareMaxEnt12(...
  maxent1_large,maxent2_large,freq_large,...
  'first-order maxent','second-order maxent','raw freqs');

% save
save(...
  sprintf('data/standard/%iingred/large.mat',Ningred_used),...
  'maxent1_large',...
  'maxent2_large',...
  'freq_large',...
  'maxent1v2',...
  'Ningred',...
  'Nrec',...
  'Nsubsamps',...
  'subsamp_size',...
  'Ningred_used',...
  'ingreds_subsamp',...
  'recipes_binary_subsamp',...
  '-v7.3');
save(...
  sprintf('data/standard/%iingred/small.mat',Ningred_used),...
  'maxent1_small',...
  'maxent2_small',...
  'freq_small',...
  'maxent1v2',...
  'Ningred',...
  'Nrec',...
  'Nsubsamps',...
  'subsamp_size',...
  'Ningred_used',...
  '-v7.3');

end