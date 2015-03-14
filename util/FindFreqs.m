function [output_small,output_large] = FindFreqs(data,subsamp_size,Nsubsamp)

% Takes a binary matrix whose cols represent data samples at
% given points in time and estimates pattern frequencies with error bars.
% Pattern frequencies and error bars are estimated by subsampling the data
% to size subsamp_size Nsubsamp times, calculating pattern
% frequencies for each subsampling, and taking the mean and standard
% deviation of each pattern frequency.
% 
% Returns structure whose entries are the mean pattern freqs, standard
% deviations of the pattern freqs, the total number of samples in the
% original data set, the size of each subsampled data set, the number of
% subsampled data sets, and a matrix containing the pattern freqs in each
% subsampled data set.
% 
% Written by: DJ Strouse
% Last updated: Aug 19, 2013 by DJ Strouse
% Part of: MaxEnt code suite
% 
% INPUTS
% data [=] data_dim X samp_size
% subsamp_size [=] positive integer
% Nsubsamp [=] integer>=1
% 
% OUTPUTS
% output_large [=] struct containing all outputs
% output_small [=] struct containing minimal outputs
% ----- if Nsubsamp==1, the above are the *same* -----
% output.probs [=] (2^data_dim) X 1
% output.samp_size [=] positive integer
% output.subsamp_size [=] positive integer
% output.Nsubsamp [=] positive integer
% output.probs_10e6 [=] (2^data_dim) X Nsubsamp
% output.probs_10e2 [=] (2^data_dim) X Nsubsamp
% output.N10e6s [=] Nsubsamp X 1
% output.N10e2s [=] Nsubsamp X 1
% output.ents [=] Nsubsamp X 1
% output.fitting_times [=] Nsubsamp X 1
% ----- if Nsubsamp>1 -----
% output.prob_ebs [=] (2^data_dim) X 1
% output.subsamp_probs [=] (2^data_dim) X Nsubsamp
% output.N10e6_mean [=] non-negative real
% output.N10e2_mean [=] non-negative real
% output.N10e6_std [=] non-negative real
% output.N10e2_std [=] non-negative real
% output.ent_mean [=] non-negative real
% output.ent_std [=] non-negative real

% init
if ~exist('Nsubsamp','var')
  disp('defaulted to 25 subsamplings')
  Nsubsamp = 25; % default to 25 subsamplings
end
if ~exist('subsamp_size','var')
  subsamp_size = round(.5*size(data,2)); % default to subsample fraction of 50%
  disp(sprintf('defaulted to subsamps of size %i (half of data)',subsamp_size))
end
data_dim = size(data,1);
samp_size = size(data,2);
num_patterns = 2^data_dim;
subsamp_probs = zeros(num_patterns, Nsubsamp);
fitting_times = zeros(Nsubsamp,1);

for i = 1:Nsubsamp
    
  tic;
  disp(sprintf('beginning subsampling %i of %i',i,Nsubsamp))

  % randomly delete subsamp_frac of the data
  subsamp_data = data(:,randperm(samp_size,subsamp_size));

  % check dims
  if size(subsamp_data,2)~=subsamp_size
    error('Dimensions of subsamp_data incompatible with subsamp_size!')
  end

  % calculate number of occurrences of each pattern in subsampled data
  for j = 1:subsamp_size
    patt_id = bin2int(subsamp_data(:,j))+1;
    subsamp_probs(patt_id,i) = subsamp_probs(patt_id,i)+1;
  end
  clear j;

  % report
  fitting_times(i) = toc;
  disp(sprintf('finished subsampling %i of %i in %.2f minutes.',...
    i,Nsubsamp,fitting_times(i)/60))

end
clear i;

% normalize
subsamp_probs = subsamp_probs/subsamp_size;

% calculate mean and standard deviation of freqs
probs = mean(subsamp_probs,2);
if Nsubsamp>1
  prob_ebs = std(subsamp_probs')';
end

% calculate entropy and identify high probability patterns
probs_10e6 = subsamp_probs>10^-6; % mask for probabilities larger than 10e-6
probs_10e2 = subsamp_probs>10^-2; % mask for probabilities larger than 10e-2
N10e6s = sum(probs_10e6,1)';
N10e2s = sum(probs_10e2,1)';
if Nsubsamp>1
  N10e6_mean = mean(N10e6s);
  N10e2_mean = mean(N10e2s);
  N10e6_std = std(N10e6s);
  N10e2_std = std(N10e2s);
end
ents = entropy(subsamp_probs);
if Nsubsamp>1
  ent_mean = mean(ents);
  ent_std = std(ents);
end

if Nsubsamp>1
  output_large = struct(...
    'probs',probs,...
    'prob_ebs',prob_ebs,...
    'samp_size',samp_size,...
    'subsamp_size',subsamp_size,...
    'Nsubsamp',Nsubsamp,...
    'subsamp_probs',subsamp_probs,...
    'probs_10e6',probs_10e6,...
    'probs_10e2',probs_10e2,...
    'N10e6s',N10e6s,...
    'N10e2s',N10e2s,...
    'N10e6_mean',N10e6_mean,...
    'N10e2_mean',N10e2_mean,...
    'N10e6_std',N10e6_std,...
    'N10e2_std',N10e2_std,...
    'ents',ents,...
    'ent_mean',ent_mean,...
    'ent_std',ent_std,...
    'fitting_times',fitting_times);
  output_small = struct(...
    'probs',probs,...
    'prob_ebs',prob_ebs,...
    'samp_size',samp_size,...
    'subsamp_size',subsamp_size,...
    'Nsubsamp',Nsubsamp,...
    'ent_mean',ent_mean,...
    'ent_std',ent_std,...
    'fitting_times',fitting_times);
else
  output_large = struct(...
    'probs',probs,...
    'samp_size',samp_size,...
    'subsamp_size',subsamp_size,...
    'Nsubsamp',Nsubsamp,...
    'probs_10e6',probs_10e6,...
    'probs_10e2',probs_10e2,...
    'N10e6s',N10e6s,...
    'N10e2s',N10e2s,...
    'ents',ents,...
    'fitting_times',fitting_times);
  output_small = output_large;
end

end