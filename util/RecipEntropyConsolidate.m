function [] = RecipEntropyConsolidate(Ningred_used,orders,Nsubsamps)
% Consolidates subsamples for order-order maxent model of Ningred_used most
% common ingredients
%
% Written by: DJ Strouse
% Last updated: Aug 19, 2013 by DJ Strouse
%
% INPUTS
% Ningred_used [=] positive integer = number of ingredients to use
% orders [=] vector of non-negative integers = orders of maxent models (0=freqs)
% Nsubsamps [=] vector of positive integers = Nsubsamps for orders
%
% OUTPUTS
% none

% init
if ~exist('orders','var')
  orders = 0:2; % defaults to cycling through all freqs and 1st/2nd maxent
end
if ~exist('Nsubsamps','var')
  Nsubsamps = [25 25 25]; % defaults to 25 subsamplings for all orders
end
if (length(Nsubsamps)==1)&&(length(orders)>1)
  % if only one Nsubsamps, copy for all orders
  Nsubsamps = Nsubsamps*ones(length(orders),1);
end
if length(orders)~=length(Nsubsamps)
  error('Lengths of orders and Nsubsamps inconsistent!')
end
Ncombos = round(2^Ningred_used);

% cycle through specified orders
for o = orders
  % change to appropriate directory
  if (length(orders)>1)&&(o~=orders(1)) % if multiple orders and after first
    cd('../../../../'); % then remember to pop up to original directory
  end
  cd(sprintf('data/standard/%iingred/order%i',Ningred_used,o));
  if o==0 % raw freqs
    % init data structures
    subsamp_probs = zeros(Ncombos,Nsubsamps(o+1));
    fitting_times = zeros(Nsubsamps(o+1),1);
    for s = 1:Nsubsamps(o+1)
      % load data
      load(sprintf('subsamp%i/subsample.mat',s));
      % extract subsample data
      subsamp_probs(:,s) = freq.probs;
      fitting_times(s) = freq.fitting_times;
      % check for consistency across subsamplings
      if s==1
        samp_size = freq.samp_size;
        subsamp_size = freq.subsamp_size;
      else
        if samp_size~=freq.samp_size
          error('samp_size inconsistent across subsamplings!')
        end
        if subsamp_size~=freq.subsamp_size
          error('subsamp_size inconsistent across subsamplings!')
        end
      end
      clear freq;
    end
    clear s;
    % calculate mean and standard deviation of freqs
    probs = mean(subsamp_probs,2);
    prob_ebs = std(subsamp_probs,0,2);
    % calculate entropy and identify high probability patterns
    probs_10e6 = subsamp_probs>10^-6; % mask for probabilities larger than 10e-6
    probs_10e2 = subsamp_probs>10^-2; % mask for probabilities larger than 10e-2
    N10e6s = sum(probs_10e6,1)';
    N10e2s = sum(probs_10e2,1)';
    N10e6_mean = mean(N10e6s);
    N10e2_mean = mean(N10e2s);
    N10e6_std = std(N10e6s);
    N10e2_std = std(N10e2s);
    ents = entropy(subsamp_probs);
    ent_mean = mean(ents);
    ent_std = std(ents);
    % package and save
    freq_large = struct(...
      'probs',probs,...
      'prob_ebs',prob_ebs,...
      'samp_size',samp_size,...
      'subsamp_size',subsamp_size,...
      'Nsubsamp',Nsubsamps(o+1),...
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
    freq_small = struct(...
      'probs',probs,...
      'prob_ebs',prob_ebs,...
      'samp_size',samp_size,...
      'subsamp_size',subsamp_size,...
      'Nsubsamp',Nsubsamps(o+1),...
      'ent_mean',ent_mean,...
      'ent_std',ent_std,...
      'fitting_times',fitting_times);
    clear probs prob_ebs samp_size subsamp_size subsamp_probs;
    clear probs_10e6 probs_10e2 N10e6s N10e2s N10e6s_mean N10e2s_mean;
    clear N10e6s_std N10e2s_std ents ents_mean ents_std fitting_times;
    save('large.mat','freq_large','-v7.3');
    save('small.mat','freq_small','-v7.3');
    clear freq_large freq_small;
  else % maxent model
    % init data structures
    subsamp_probs = zeros(Ncombos,Nsubsamps(o+1));
    subsamp_param = cell(o,1);
    for i = 1:o
      subsamp_param{i} = zeros(nchoosek(Ningred_used,i),Nsubsamps(o+1));
    end
    clear i;
    subsamp_Z = ones(Nsubsamps(o+1),1);
    fitting_times = zeros(Nsubsamps(o+1),1);
    for s = 1:Nsubsamps(o+1)
      % load data
      load(sprintf('subsamp%i/subsample.mat',s));
      % extract subsample data
      subsamp_probs(:,s) = maxent.probs;
      fitting_times(s) = maxent.fitting_times;
      subsamp_Z(s) = maxent.Z;
      if o>1 % fixes problem with structs converting single-el cells to arrays
        for i = 1:o
          subsamp_param{i}(:,s) = maxent(i).param;
        end
        clear i;
      else
        subsamp_param{1}(:,s) = maxent.param;
      end
      % check for consistency across subsamplings
      if s==1
        samp_size = maxent(1).samp_size;
        subsamp_size = maxent(1).subsamp_size;
      else
        if samp_size~=maxent(1).samp_size
          error('samp_size inconsistent across subsamplings!')
        end
        if subsamp_size~=maxent(1).subsamp_size
          error('subsamp_size inconsistent across subsamplings!')
        end
      end
      clear maxent;
    end
    clear s;
    % calculate mean/std of estimated parameters for output
    probs = mean(subsamp_probs,2);
    param = cell(o,1);
    prob_ebs = std(subsamp_probs,0,2);
    param_ebs = cell(o,1);
    for i = 1:o
      param{i} = mean(subsamp_param{i},2);
      param_ebs{i} = std(subsamp_param{i},0,2);
    end
    clear i;
    Z = mean(subsamp_Z,2);
    Z_std = std(subsamp_Z,0,2);

    % calculate entropy and identify high probability patterns
    probs_10e6 = subsamp_probs>10^-6; % mask for probabilities larger than 10e-6
    probs_10e2 = subsamp_probs>10^-2; % mask for probabilities larger than 10e-2
    N10e6s = sum(probs_10e6,1)';
    N10e2s = sum(probs_10e2,1)';
    N10e6_mean = mean(N10e6s);
    N10e2_mean = mean(N10e2s);
    N10e6_std = std(N10e6s);
    N10e2_std = std(N10e2s);
    ents = entropy(subsamp_probs);
    ent_mean = mean(ents);
    ent_std = std(ents);

    % unpack param mean and ebs for more intuitive handling (vector -> arrays)
    param_array = cell(o,1); % each entry is an array for order index
    param_ebs_array = cell(o,1);
    indices = GenSubsetIndices(Ningred_used,o);
    for i = 1:o % this is about to get hacky as shit - make pretty later
      if i==1 % copy parameters for first-order
        param_array{i} = param{i}; % Ningred_used X 1
        param_ebs_array{i} = param_ebs{i}; % Ningred_used X 1
      else
        if i==2 % copy parameters for second-order
          param_array{i} = zeros(Ningred_used,Ningred_used);
          param_ebs_array{i} = zeros(Ningred_used,Ningred_used);
          for j = 1:length(param{i})
            param_array{i}(indices{i}(1,j),indices{i}(2,j)) = param{i}(j);
            param_ebs_array{i}(indices{i}(1,j),indices{i}(2,j)) = param_ebs{i}(j);
          end
          clear j;
        elseif i==3 % copy parameters for third-order
          param_array{i} = zeros(Ningred_used,Ningred_used,Ningred_used);
          param_ebs_array{i} = zeros(Ningred_used,Ningred_used,Ningred_used);
          for j = 1:length(param{i})
            param_array{i}(...
              indices{i}(1,j),indices{i}(2,j),indices{i}(3,j)...
              ) = param{i}(j);
            param_ebs_array{i}(...
              indices{i}(1,j),indices{i}(2,j),indices{i}(3,j)...
              ) = param_ebs{i}(j);
          end
          clear j;
        elseif i==4 % copy parameters for fourth-order
          param_array{i} = zeros(Ningred_used,Ningred_used,Ningred_used,Ningred_used);
          param_ebs_array{i} = zeros(Ningred_used,Ningred_used,Ningred_used,Ningred_used);
          for j = 1:length(param{i})
            param_array{i}(...
              indices{i}(1,j),indices{i}(2,j),indices{i}(3,j),indices{i}(4,j)...
              ) = param{i}(j);
            param_ebs_array{i}(...
              indices{i}(1,j),indices{i}(2,j),indices{i}(3,j),indices{i}(4,j)...
              ) = param_ebs{i}(j);
          end
          clear j;
        end
        param_array{i} = symmetrize(param_array{i});
        param_ebs_array{i} = symmetrize(param_ebs_array{i});
      end
    end
    clear i;

    % package and save
    maxent_large = struct(...
      'order',o,...
      'probs',probs,...
      'prob_ebs',prob_ebs,...
      'samp_size',samp_size,...
      'subsamp_size',subsamp_size,...
      'Nsubsamp',Nsubsamps(o+1),...
      'param',param,...
      'param_ebs',param_ebs,...
      'param_array',param_array,...
      'param_ebs_array',param_ebs_array,...
      'Z',Z,...
      'Z_std',Z_std,...
      'subsamp_probs',subsamp_probs,...
      'subsamp_param',subsamp_param,...
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
    maxent_small = struct(...
      'order',o,...
      'probs',probs,...
      'prob_ebs',prob_ebs,...
      'samp_size',samp_size,...
      'subsamp_size',subsamp_size,...
      'Nsubsamp',Nsubsamps(o+1),...
      'ent_mean',ent_mean,...
      'ent_std',ent_std,...
      'fitting_times',fitting_times);
    clear probs prob_ebs samp_size subsamp_size subsamp_probs;
    clear probs_10e6 probs_10e2 N10e6s N10e2s N10e6s_mean N10e2s_mean;
    clear N10e6s_std N10e2s_std ents ents_mean ents_std fitting_times;
    clear param param_ebs param_array param_ebs_array Z Z_std subsamp_param;
    save('large.mat','maxent_large','-v7.3');
    save('small.mat','maxent_small','-v7.3');
    clear maxent_large maxent_small;
  end
end
clear o;
cd('../../../../');

end

