function [] = RecipEntropySpicesSingle(Ningred_used,order,subsamp_index,save_flag)
% Fits a single maxent model of order order using the Ningred_used most
% common spices and saves it as subsample number subsamp_index
%
% Written by: DJ Strouse
% Last updated: Aug 22, 2013 by DJ Strouse
%
% INPUTS
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
load('data/spices.mat');
recipes_binary = final_spices_binary; clear final_spices_binary;
ingreds = final_spices; clear final_spices;
mkdir(sprintf('data/spices/%iingred/order%i/subsamp%i',...
  Ningred_used,order,subsamp_index)); % creates appropriate directory
cd(sprintf('data/spices/%iingred/order%i/subsamp%i',...
  Ningred_used,order,subsamp_index)); % changes to appropriate directory
diary('diary.txt');
Ningred = size(recipes_binary,2);
Nrec = size(recipes_binary,1);
subsamp_size = round(.8*Nrec);
ingred_freq = ingred_freq_spices;
if Ningred_used>Ningred
  error('Ningred_used>Ningred!')
end
disp(sprintf('Subsampling from %i recipes down to %i',Nrec,subsamp_size))

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