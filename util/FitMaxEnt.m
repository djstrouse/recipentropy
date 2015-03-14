function [output_small,output_large] =...
  FitMaxEnt(data,order,subsamp_size,Nsubsamp,method,progress)

% Takes a binary matrix whose cols represent data samples at
% given points in time and fits a MaxEnt model to the pattern probabilities
% using correlations up to order. Error bars are also estimated by fitting
% multiple MaxEnt models to subsamplings of the data (and using the mean
% and standard deviation of the estimated parameters).
% May periodically save progress depending on single subsampling fitting
% times. Takes an optional argument which specifies fitting method: exact
% or minimum probability flow (MPF). The latter should be used for any data
% sets with ~20+ variables (but only works for order one and two MaxEnt
% models. Also takes an optional argument which contains progress from a
% previous fitting and, if input, contains from that point.
% 
%
% UPDATE THIS FOR MPF
% Returns structure of probabilities (and error bars), number of samples
% in original data set, size of subsampled data sets, the number of
% subsampled data sets, order parameters (and error bars) in vector and
% array form, normalization constant Z (and error bars), matrices
% containing parameters and probabilities estimated on each subsampled
% data set, and other statistics (e.g. entropy).
% 
% Note: currently works only up to order four. Can easily be tidied up to
% accommodate higher-order maxent models in the future.
% 
% Written by: DJ Strouse
% Last updated: Feb 8, 2014 by DJ Strouse
% Part of: MaxEnt code suite
% 
% INPUTS
% data [=] data_dim X samp_size
% order [=] positive integer (1 or 2 if method='mpf', 1-4 for 'exact')
% subsamp_size [=] positive integer
% Nsubsamp [=] positive integer
% method [=] string in {'exact','mpf'}
% progress [=] struct containing progress of partially finished fitting
% 
% OUTPUTS
% output_large [=] struct containing all outputs
% output_small [=] struct containing minimal outputs
% ----- if Nsubsamp==1, the above are the *same* -----
% output.order [=] positive integer
% output.probs [=] (2^data_dim) X 1
% output.samp_size [=] positive integer
% output.subsamp_size [=] positive integer
% output.Nsubsamp [=] positive integer
% output.param [=] order X 1 cell array
% output.param{j} [=] (data_dim choose j) X 1
% output.param_array [=] order X 1 cell array
% output.param_array{j} [=] data_dim^j (e.g. if j=2, data_dim X data_dim)
% output.Z [=] positive real
% output.probs_10e6 [=] (2^data_dim) X Nsubsamp
% output.probs_10e2 [=] (2^data_dim) X Nsubsamp
% output.N10e6s [=] Nsubsamp X 1
% output.N10e2s [=] Nsubsamp X 1
% output.ents [=] Nsubsamp X 1
% output.fitting_times [=] Nsubsamp X 1
% ----- if Nsubsamp>1 -----
% output.prob_ebs [=] (2^data_dim) X 1
% output.param_ebs [=] order X 1 cell array
% output.param_ebs{j} [=] (data_dim choose j) X 1
% output.param_ebs_array [=] order X 1 cell array
% output.param_ebs_array{j} [=] data_dim^j
% output.Z_std [=] positive real
% output.subsamp_probs [=] (2^data_dim) X Nsubsamp
% output.subsamp_param [=] order X 1 cell array
% output.subsamp_param{j} [=] (data_dim choose j) X Nsubsamp
% output.N10e6_mean [=] non-negative real
% output.N10e2_mean [=] non-negative real
% output.N10e6_std [=] non-negative real
% output.N10e2_std [=] non-negative real
% output.ent_mean [=] non-negative real
% output.ent_std [=] non-negative real
% 
% FUTURE IMPROVEMENTS
% better initialization
% Hessian-aided optimization
% user control over periodic saving

disp(sprintf('beginning order %i maxent fitting',order))

% read method or default to exact
if ~exist('method','var')
  method = 'exact';
elseif strcmp(method,'exact')
  if (order~=1)&&(order~=2)&&(order~=3)&&(order~=4)
    error('order must be integer 1-4 for "exact" method!')
  end
elseif strcmp(method,'exact')
  if (order~=1)&&(order~=2)
    error('order must be integer 1-2 for "mpf" method!')
  end
end

% exact fitting method
if strcmp(method,'exact')
  
  % check for progress input
  if exist('progress','var')

    first_subsamp = progress.last_subsamp;
    order = progress.order;
    samp_size = progress.samp_size;
    subsamp_size = progress.subsamp_size;
    Nsubsamp = progress.Nsubsamp;
    subsamp_probs = progress.subsamp_probs;
    subsamp_param = progress.subsamp_param;
    fitting_times = progress.fitting_times;

  else

    % init
    first_subsamp = 1;
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
    fitting_times = zeros(Nsubsamp,1);

    % calculate number of parameters required
    num_param = data_dim;
    for i = 2:order
      num_param = num_param+nchoosek(data_dim,i);
    end
    clear i;

    % initializations that will not change during subsampling
    lambda0 = ones(num_param,1); % may want to toy with different initial guesses
    [indices,subsets] = GenSubsetIndices(data_dim,order);
    f_data = f(data,order,subsets); % pre-processes data
    options = optimset('GradObj','on'); % helpful to calc/include a Hessian?
    subsamp_param = cell(order,1);
    for i = 1:order
      subsamp_param{i} = zeros(nchoosek(data_dim,i),Nsubsamp);
    end
    clear i;
    subsamp_probs = ones(num_patterns,Nsubsamp);
    subsamp_Z = ones(Nsubsamp,1);

  end

  for i = first_subsamp:Nsubsamp

    tic;
    disp(sprintf('beginning subsampling %i of %i',i,Nsubsamp))

    % randomly delete subsamp_frac of the data in f_data
    subsamp_f_data = f_data(:,randperm(samp_size,subsamp_size));

    % check dims
    if size(subsamp_f_data,2)~=subsamp_size
      error('Dimensions of subsamp_f_data incompatible with subsamp_size!')
    end

    % initialize parameters for max likelihood
    data_avgd_subsamp_f = (1/subsamp_size)*sum(subsamp_f_data,2);
                                       % avgs subsamp_f_data across samples
    obj_func = @(lambda)(NegLogLik(data_avgd_subsamp_f,subsamp_size,...
        data_dim,order,subsets,lambda));

    % build MaxEnt model on subsampled data
    lambda_hat = fminunc(obj_func,lambda0,options);

    % store parameters and probabilities

    % separately store parameters for different orders
    % cell j indicates order j parameters
    % within a cell, columns correspond to different subsamplings while
    % rows correspond to different parameters
    count = 1;
    for j = 1:order
      num_param = nchoosek(data_dim,j);
      subsamp_param{j}(:,i) = lambda_hat(count:num_param+count-1);
      count = count+num_param;
    end
    clear j;

    % calculate unnormalized probabilities
    % columns correspond to different subsamplings while rows correspond to
    % different probabilities
    for j = 1:num_patterns
      pattern = int2bin(j-1,data_dim);
      f_pattern = f(pattern',order,subsets);
      subsamp_probs(j,i) = exp(f_pattern'*lambda_hat);
    end
    clear j;

    % calculate norm const Z
    subsamp_Z(i) = sum(subsamp_probs(:,i));

    % normalize probabilities
    subsamp_probs(:,i) = subsamp_probs(:,i)/subsamp_Z(i);

    % report
    fitting_times(i) = toc;
    disp(sprintf('finished subsampling %i of %i in %.2f minutes.',...
      i,Nsubsamp,fitting_times(i)/60))

    % save
    if ((fitting_times(i)>(2*60*60)&&...           % if fitting time >2 hours
        mod(i,5)==0)||...              %...and subsampling is a multiple of 5
        fitting_times(i)>(12*60*60))&&...    %...or if fitting time >12 hours
        i~=Nsubsamp                  %...and this is not the last subsampling
      progress = struct(...                            %...then save progress
        'last_subsamp',i,...
        'order',order,...
        'samp_size',samp_size,...
        'subsamp_size',subsamp_size,...
        'Nsubsamp',Nsubsamp,...
        'subsamp_probs',subsamp_probs,...
        'subsamp_param',subsamp_param,...
        'fitting_times',fitting_times);
      save(...
        sprintf('maxent%i_subsampling%iof%i',order,i,Nsubsamp),...
        'progress',...
        '-v7.3');
      disp('successfully saved progress')
    end

  end
  clear i;

  % calculate mean/std of estimated parameters for output
  probs = mean(subsamp_probs,2);
  param = cell(order,1);
  if Nsubsamp>1
    prob_ebs = std(subsamp_probs')';
    param_ebs = cell(order,1);
  end
  for i = 1:order
    param{i} = mean(subsamp_param{i},2);
    if Nsubsamp>1
      param_ebs{i} = std(subsamp_param{i}')';
    end
  end
  clear i;
  Z = mean(subsamp_Z,2);
  if Nsubsamp>1
    Z_std = std(subsamp_Z,0,2);
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

  % unpack param mean and ebs for more intuitive handling (vector -> arrays)
  param_array = cell(order,1); % each entry is an array for order index
  if Nsubsamp>1
    param_ebs_array = cell(order,1);
  end
  for i = 1:order % this is about to get hacky as shit - make pretty later
    if i==1 % copy parameters for first-order
      param_array{i} = param{i}; % data_dim X 1
      if Nsubsamp>1
        param_ebs_array{i} = param_ebs{i}; % data_dim X 1
      end
    else
      if i==2 % copy parameters for second-order
        param_array{i} = zeros(data_dim,data_dim);
        if Nsubsamp>1
          param_ebs_array{i} = zeros(data_dim,data_dim);
        end
        for j = 1:length(param{i})
          param_array{i}(indices{i}(1,j),indices{i}(2,j)) = param{i}(j);
          if Nsubsamp>1
            param_ebs_array{i}(indices{i}(1,j),indices{i}(2,j)) = param_ebs{i}(j);
          end
        end
        clear j;
      elseif i==3 % copy parameters for third-order
        param_array{i} = zeros(data_dim,data_dim,data_dim);
        if Nsubsamp>1
          param_ebs_array{i} = zeros(data_dim,data_dim,data_dim);
        end
        for j = 1:length(param{i})
          param_array{i}(...
            indices{i}(1,j),indices{i}(2,j),indices{i}(3,j)...
            ) = param{i}(j);
          if Nsubsamp>1
            param_ebs_array{i}(...
              indices{i}(1,j),indices{i}(2,j),indices{i}(3,j)...
              ) = param_ebs{i}(j);
          end
        end
        clear j;
      elseif i==4 % copy parameters for fourth-order
        param_array{i} = zeros(data_dim,data_dim,data_dim,data_dim);
        if Nsubsamp>1
          param_ebs_array{i} = zeros(data_dim,data_dim,data_dim,data_dim);
        end
        for j = 1:length(param{i})
          param_array{i}(...
            indices{i}(1,j),indices{i}(2,j),indices{i}(3,j),indices{i}(4,j)...
            ) = param{i}(j);
          if Nsubsamp>1
            param_ebs_array{i}(...
              indices{i}(1,j),indices{i}(2,j),indices{i}(3,j),indices{i}(4,j)...
              ) = param_ebs{i}(j);
          end
        end
        clear j;
      end
      param_array{i} = symmetrize(param_array{i});
      if Nsubsamp>1
        param_ebs_array{i} = symmetrize(param_ebs_array{i});
      end
    end
    clear i;

    if Nsubsamp>1
      output_large = struct(...
        'order',order,...
        'probs',probs,...
        'prob_ebs',prob_ebs,...
        'samp_size',samp_size,...
        'subsamp_size',subsamp_size,...
        'Nsubsamp',Nsubsamp,...
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
      output_small = struct(...
        'order',order,...
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
        'order',order,...
        'probs',probs,...
        'samp_size',samp_size,...
        'subsamp_size',subsamp_size,...
        'Nsubsamp',Nsubsamp,...
        'param',param,...
        'param_array',param_array,...
        'Z',Z,...
        'probs_10e6',probs_10e6,...
        'probs_10e2',probs_10e2,...
        'N10e6s',N10e6s,...
        'N10e2s',N10e2s,...
        'ents',ents,...
        'fitting_times',fitting_times);
      output_small = output_large;
    end

    disp(sprintf('finished order %i maxent fitting',order))

  end

% mpf fitting method
elseif strcmp(method,'mpf')
  
  % check for progress input
  if exist('progress','var')

    % COME BACK AND FIX THIS
%     first_subsamp = progress.last_subsamp;
%     order = progress.order;
%     samp_size = progress.samp_size;
%     subsamp_size = progress.subsamp_size;
%     Nsubsamp = progress.Nsubsamp;
%     subsamp_probs = progress.subsamp_probs;
%     subsamp_param = progress.subsamp_param;
%     fitting_times = progress.fitting_times;

  else

    % init
    first_subsamp = 1;
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
    fitting_times = zeros(Nsubsamp,1);

    % calculate number of parameters required
    num_param = data_dim;
    for i = 2:order
      num_param = num_param+nchoosek(data_dim,i);
    end
    clear i;
    % HERE HERE HERE HERE HERE HERE HERE HERE HERE HERE
    % initializations that will not change during subsampling
    lambda0 = ones(num_param,1); % may want to toy with different initial guesses
    [indices,subsets] = GenSubsetIndices(data_dim,order);
    f_data = f(data,order,subsets); % pre-processes data
    options = optimset('GradObj','on'); % helpful to calc/include a Hessian?
    subsamp_param = cell(order,1);
    for i = 1:order
      subsamp_param{i} = zeros(nchoosek(data_dim,i),Nsubsamp);
    end
    clear i;
    subsamp_probs = ones(num_patterns,Nsubsamp);
    subsamp_Z = ones(Nsubsamp,1);

  end

  for i = first_subsamp:Nsubsamp

    tic;
    disp(sprintf('beginning subsampling %i of %i',i,Nsubsamp))

    % randomly delete subsamp_frac of the data in f_data
    subsamp_f_data = f_data(:,randperm(samp_size,subsamp_size));

    % check dims
    if size(subsamp_f_data,2)~=subsamp_size
      error('Dimensions of subsamp_f_data incompatible with subsamp_size!')
    end

    % initialize parameters for max likelihood
    data_avgd_subsamp_f = (1/subsamp_size)*sum(subsamp_f_data,2);
                                       % avgs subsamp_f_data across samples
    obj_func = @(lambda)(NegLogLik(data_avgd_subsamp_f,subsamp_size,...
        data_dim,order,subsets,lambda));

    % build MaxEnt model on subsampled data
    lambda_hat = fminunc(obj_func,lambda0,options);

    % store parameters and probabilities

    % separately store parameters for different orders
    % cell j indicates order j parameters
    % within a cell, columns correspond to different subsamplings while
    % rows correspond to different parameters
    count = 1;
    for j = 1:order
      num_param = nchoosek(data_dim,j);
      subsamp_param{j}(:,i) = lambda_hat(count:num_param+count-1);
      count = count+num_param;
    end
    clear j;

    % calculate unnormalized probabilities
    % columns correspond to different subsamplings while rows correspond to
    % different probabilities
    for j = 1:num_patterns
      pattern = int2bin(j-1,data_dim);
      f_pattern = f(pattern',order,subsets);
      subsamp_probs(j,i) = exp(f_pattern'*lambda_hat);
    end
    clear j;

    % calculate norm const Z
    subsamp_Z(i) = sum(subsamp_probs(:,i));

    % normalize probabilities
    subsamp_probs(:,i) = subsamp_probs(:,i)/subsamp_Z(i);

    % report
    fitting_times(i) = toc;
    disp(sprintf('finished subsampling %i of %i in %.2f minutes.',...
      i,Nsubsamp,fitting_times(i)/60))

    % save
    if ((fitting_times(i)>(2*60*60)&&...           % if fitting time >2 hours
        mod(i,5)==0)||...              %...and subsampling is a multiple of 5
        fitting_times(i)>(12*60*60))&&...    %...or if fitting time >12 hours
        i~=Nsubsamp                  %...and this is not the last subsampling
      progress = struct(...                            %...then save progress
        'last_subsamp',i,...
        'order',order,...
        'samp_size',samp_size,...
        'subsamp_size',subsamp_size,...
        'Nsubsamp',Nsubsamp,...
        'subsamp_probs',subsamp_probs,...
        'subsamp_param',subsamp_param,...
        'fitting_times',fitting_times);
      save(...
        sprintf('maxent%i_subsampling%iof%i',order,i,Nsubsamp),...
        'progress',...
        '-v7.3');
      disp('successfully saved progress')
    end

  end
  clear i;

  % calculate mean/std of estimated parameters for output
  probs = mean(subsamp_probs,2);
  param = cell(order,1);
  if Nsubsamp>1
    prob_ebs = std(subsamp_probs')';
    param_ebs = cell(order,1);
  end
  for i = 1:order
    param{i} = mean(subsamp_param{i},2);
    if Nsubsamp>1
      param_ebs{i} = std(subsamp_param{i}')';
    end
  end
  clear i;
  Z = mean(subsamp_Z,2);
  if Nsubsamp>1
    Z_std = std(subsamp_Z,0,2);
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

  % unpack param mean and ebs for more intuitive handling (vector -> arrays)
  param_array = cell(order,1); % each entry is an array for order index
  if Nsubsamp>1
    param_ebs_array = cell(order,1);
  end
  for i = 1:order % this is about to get hacky as shit - make pretty later
    if i==1 % copy parameters for first-order
      param_array{i} = param{i}; % data_dim X 1
      if Nsubsamp>1
        param_ebs_array{i} = param_ebs{i}; % data_dim X 1
      end
    else
      if i==2 % copy parameters for second-order
        param_array{i} = zeros(data_dim,data_dim);
        if Nsubsamp>1
          param_ebs_array{i} = zeros(data_dim,data_dim);
        end
        for j = 1:length(param{i})
          param_array{i}(indices{i}(1,j),indices{i}(2,j)) = param{i}(j);
          if Nsubsamp>1
            param_ebs_array{i}(indices{i}(1,j),indices{i}(2,j)) = param_ebs{i}(j);
          end
        end
        clear j;
      elseif i==3 % copy parameters for third-order
        param_array{i} = zeros(data_dim,data_dim,data_dim);
        if Nsubsamp>1
          param_ebs_array{i} = zeros(data_dim,data_dim,data_dim);
        end
        for j = 1:length(param{i})
          param_array{i}(...
            indices{i}(1,j),indices{i}(2,j),indices{i}(3,j)...
            ) = param{i}(j);
          if Nsubsamp>1
            param_ebs_array{i}(...
              indices{i}(1,j),indices{i}(2,j),indices{i}(3,j)...
              ) = param_ebs{i}(j);
          end
        end
        clear j;
      elseif i==4 % copy parameters for fourth-order
        param_array{i} = zeros(data_dim,data_dim,data_dim,data_dim);
        if Nsubsamp>1
          param_ebs_array{i} = zeros(data_dim,data_dim,data_dim,data_dim);
        end
        for j = 1:length(param{i})
          param_array{i}(...
            indices{i}(1,j),indices{i}(2,j),indices{i}(3,j),indices{i}(4,j)...
            ) = param{i}(j);
          if Nsubsamp>1
            param_ebs_array{i}(...
              indices{i}(1,j),indices{i}(2,j),indices{i}(3,j),indices{i}(4,j)...
              ) = param_ebs{i}(j);
          end
        end
        clear j;
      end
      param_array{i} = symmetrize(param_array{i});
      if Nsubsamp>1
        param_ebs_array{i} = symmetrize(param_ebs_array{i});
      end
    end
    clear i;

    if Nsubsamp>1
      output_large = struct(...
        'order',order,...
        'probs',probs,...
        'prob_ebs',prob_ebs,...
        'samp_size',samp_size,...
        'subsamp_size',subsamp_size,...
        'Nsubsamp',Nsubsamp,...
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
      output_small = struct(...
        'order',order,...
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
        'order',order,...
        'probs',probs,...
        'samp_size',samp_size,...
        'subsamp_size',subsamp_size,...
        'Nsubsamp',Nsubsamp,...
        'param',param,...
        'param_array',param_array,...
        'Z',Z,...
        'probs_10e6',probs_10e6,...
        'probs_10e2',probs_10e2,...
        'N10e6s',N10e6s,...
        'N10e2s',N10e2s,...
        'ents',ents,...
        'fitting_times',fitting_times);
      output_small = output_large;
    end

    disp(sprintf('finished order %i maxent fitting',order))

  end

else
  error('method needs to be either "exact" or "mpf"!')
end