function [] = RecipEntropyCulturesSingle(culture,Ningred_used,order,subsamp_index,save_flag)
% Fits a single maxent model of order order using the Ningred_used most
% common ingredients from a specified culture and saves it as
% subsample number subsamp_index
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
% order [=] non-negative integer = order of maxent model (0=freqs)
% subsamp_index [=] positive integer = index of subsample
% save_flag [=] boolean = indicates whether or not to save data (def=true)
%
% OUTPUTS
% none

% init
if ~exist('save_flag','var')
  save_flag = true;
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
load('data/recipes.mat');
% creates appropriate directory
mkdir(['data/cultures/',culture_str,...
  sprintf('/%iingred/order%i/subsamp%i',Ningred_used,order,subsamp_index)]);
% changes to appropriate directory
cd(['data/cultures/',culture_str,...
  sprintf('/%iingred/order%i/subsamp%i',Ningred_used,order,subsamp_index)]);
diary('diary.txt');
Ningred = size(recipes_binary,2);
Nrec = size(recipes_binary,1);
if Ningred_used>Ningred
  error('Ningred_used>Ningred!')
end

% use recipes only of specified culture
[~,recipes_binary,Nrec_used]=...
  ExtractCulturalRecipes(culture,recipes,recipes_binary);
subsamp_size = round(.8*Nrec_used);

disp(sprintf('Subsampling from %i recipes down to %i',Nrec_used,subsamp_size))

% use Ningred_used most common ingredients only
recipes_binary_subsamp = ExtractIngred(Ningred_used,ingred_freq,recipes_binary,ingreds);

% fit model
if order==0
  disp('Subsampling raw frequencies')
  freq =...
    FindFreqs(recipes_binary_subsamp,subsamp_size,1);
  disp('Subsampled raw frequencies')
  if save_flag
    save(...
      'subsample.mat',...
      'freq',...
      '-v7.3');
  end
else
  disp(sprintf('Building model of order %i',order))
  maxent =...
    FitMaxEnt(recipes_binary_subsamp,order,subsamp_size,1);
  disp(sprintf('Built model of order %i',order))
  if save_flag
    save(...
      'subsample.mat',...
      'maxent',...
      '-v7.3');
  end
end

%% fixes a problem when running this in cambridge
! rm -f nohup.out

end