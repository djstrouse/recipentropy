function [output] = CompareMaxEnt12(model1,model2,freqs,...
  name1,name2,name0)

% Compares a 1st order MaxEnt model, 2nd order MaxEnt model, and observed
% pattern frequencies in data. Calculates and outputs entropies,
% connected information (of order 2), multi-information, and
% ratios of these quantities.
% 
% model1, model2, and freqs are structures containing (at least) the
% following fields:
% 'subsamp_probs' = matrix of probabilities estimated on data subsamplings
% 
% name1, name2, and name0 should be strings indicating the names of each
% model (to be displayed in output)
% 
% Written by: DJ Strouse
% Last updated: May 17, 2013 by DJ Strouse
% Part of: MaxEnt code suite
% 
% INPUTS
% model1 [=] struct (as output by fit_MaxEnt.m)
% model2 [=] struct (as output by fit_MaxEnt.m)
% freqs [=] struct (as output by find_freqs.m)
% name1 [=] string
% name2 [=] string
% name0 [=] string
% 
% OUTPUTS
% mean_H1 [=] non-negative real
% mean_H2 [=] non-negative real
% mean_H0 [=] non-negative real
% mean_I2 [=] non-negative real
% mean_IN [=] non-negative real
% mean_I2_over_H0 [=] non-negative real
% mean_IN_over_H0 [=] non-negative real
% mean_I2_over_IN [=] non-negative real
% mean_H1_over_H0 [=] non-negative real
% mean_H2_over_H0 [=] non-negative real
% std_H1 [=] non-negative real
% std_H2 [=] non-negative real
% std_H0 [=] non-negative real
% std_I2 [=] non-negative real
% std_IN [=] non-negative real
% std_I2_over_H0 [=] non-negative real
% std_IN_over_H0 [=] non-negative real
% std_I2_over_IN [=] non-negative real
% std_H1_over_H0 [=] non-negative real
% std_H2_over_H0 [=] non-negative real
% perc_error_H1 [=] non-negative real
% perc_error_H2 [=] non-negative real
% perc_error_I2 [=] non-negative real
% perc_error_IN [=] non-negative real
% perc_error_H0 [=] non-negative real
% perc_error_I2_over_H0 [=] non-negative real
% perc_error_IN_over_H0 [=] non-negative real
% perc_error_I2_over_IN [=] non-negative real
% perc_error_H1_over_H0 [=] non-negative real
% perc_error_H2_over_H0 [=] non-negative real

% calculate entropies and other statistics across subsamplings
H1 = entropy(model1(1).subsamp_probs);
H2 = entropy(model2(1).subsamp_probs);
I2 = rectify(H1-H2); % second order connected information
H0 = entropy(freqs(1).subsamp_probs);
IN = rectify(H1-H0); % multi-information
I2_over_H0 = I2./H0;
IN_over_H0 = IN./H0;
I2_over_IN = min(I2./IN,1);
H1_over_H0 = H1./H0;
H2_over_H0 = H2./H0;
for j = 1:length(I2_over_IN);
  if isnan(I2_over_IN(j))
    I2_over_IN(j) = 0;
  end
  if isinf(I2_over_IN(j))
    I2_over_IN(j) = 0;
  end
end
clear j;

% means, stds, and percent errors
mean_H1 = mean(H1);
mean_H2 = mean(H2);
mean_I2 = mean(I2);
mean_IN = mean(IN);
mean_H0 = mean(H0);
mean_I2_over_H0 = mean(I2_over_H0);
mean_IN_over_H0 = mean(IN_over_H0);
mean_I2_over_IN = mean(I2_over_IN);
mean_H1_over_H0 = mean(H1_over_H0);
mean_H2_over_H0 = mean(H2_over_H0);
std_H1 = std(H1);
std_H2 = std(H2);
std_I2 = std(I2);
std_H0 = std(H0);
std_IN = std(IN);
std_I2_over_H0 = std(I2_over_H0);
std_IN_over_H0 = std(IN_over_H0);
std_I2_over_IN = std(I2_over_IN);
std_H1_over_H0 = std(H1_over_H0);
std_H2_over_H0 = std(H2_over_H0);
perc_error_H1 = std_H1/mean_H1;
perc_error_H2 = std_H2/mean_H2;
perc_error_I2 = std_I2/mean_I2;
perc_error_IN = std_IN/mean_IN;
perc_error_H0 = std_H0/mean_H0;
perc_error_I2_over_H0 = std_I2_over_H0/mean_I2_over_H0;
perc_error_IN_over_H0 = std_IN_over_H0/mean_IN_over_H0;
perc_error_I2_over_IN = std_I2_over_IN/mean_I2_over_IN;
perc_error_H1_over_H0 = std_H1_over_H0/mean_H1_over_H0;
perc_error_H2_over_H0 = std_H2_over_H0/mean_H2_over_H0;

% display results
if std_H1~=0
  [decplace_H1 str_H1] = sigdig(std_H1); % finds n: 10^n = 1st sigdig
else % if zero std, just use all digits
  [decplace_H1 str_H1] = sigdig(mean_H1);
end
if std_H2~=0
  [decplace_H2 str_H2] = sigdig(std_H2);
else 
  [decplace_H2 str_H2] = sigdig(mean_H2);
end
if std_H0~=0
  [decplace_H0 str_H0] = sigdig(std_H0);
else
  [decplace_H0 str_H0] = sigdig(mean_H0);
end
if std_I2~=0
  [decplace_I2 str_I2] = sigdig(std_I2);
else
  [decplace_I2 str_I2] = sigdig(mean_I2);
end
if std_IN~=0
  [decplace_IN str_IN] = sigdig(std_IN);
else
  [decplace_IN str_IN] = sigdig(mean_IN);
end
if std_I2_over_H0~=0
  [decplace_I2overH0 str_I2overH0] = sigdig(std_I2_over_H0);
else
  [decplace_I2overH0 str_I2overH0] = sigdig(mean_I2_over_H0);
end
if std_IN_over_H0~=0
  [decplace_INoverH0 str_INoverH0] = sigdig(std_IN_over_H0);
else
  [decplace_INoverH0 str_INoverH0] = sigdig(mean_IN_over_H0);
end
if std_I2_over_IN~=0
  [decplace_I2overIN str_I2overIN] = sigdig(std_I2_over_IN);
else
  [decplace_I2overIN str_I2overIN] = sigdig(mean_I2_over_IN);
end
if std_H1_over_H0~=0
  [decplace_H1overH0 str_H1overH0] = sigdig(std_H1_over_H0);
else
  [decplace_H1overH0 str_H1overH0] = sigdig(mean_H1_over_H0);
end
if std_H2_over_H0~=0
  [decplace_H2overH0 str_H2overH0] = sigdig(std_H2_over_H0);
else
  [decplace_H2overH0 str_H2overH0] = sigdig(mean_H2_over_H0);
end

report = char(...
  ['entropy (H1) of ',...
  name1,...
  sprintf([' = ',str_H1],mean_H1),...
  char(177),...
  sprintf(str_H1,std_H1)],...
  ['entropy (H2) of ',...
  name2,...
  sprintf([' = ',str_H2],mean_H2),...
  char(177),...
  sprintf(str_H2,std_H2)],...
  ['entropy (H0) of ',...
  name0,...
  sprintf([' = ',str_H0],mean_H0),...
  char(177),...
  sprintf(str_H0,std_H0)],...
  [sprintf(['I2 = ',str_I2],mean_I2),...
  char(177),...
  sprintf(str_I2,std_I2)],...
  [sprintf(['IN = ',str_IN],mean_IN),...
  char(177),...
  sprintf(str_IN,std_IN)],...
  [sprintf(['I2/H0 = ',str_I2overH0],mean_I2_over_H0),...
  char(177),...
  sprintf(str_I2overH0,std_I2_over_H0)],...
  [sprintf(['IN/H0 = ',str_INoverH0],mean_IN_over_H0),...
  char(177),...
  sprintf(str_INoverH0,std_IN_over_H0)],...
  [sprintf(['I2/IN = ',str_I2overIN],mean_I2_over_IN),...
  char(177),...
  sprintf(str_I2overIN,std_I2_over_IN)],...
  [sprintf(['H1/H0 = ',str_H1overH0],mean_H1_over_H0),...
  char(177),...
  sprintf(str_H1overH0,std_H1_over_H0)],...
  [sprintf(['H2/H0 = ',str_H2overH0],mean_H2_over_H0),...
  char(177),...
  sprintf(str_H2overH0,std_H2_over_H0)]);
disp(report);

output = struct(...
  'mean_H1',mean_H1,...
  'mean_H2',mean_H2,...
  'mean_H0',mean_H0,...
  'mean_I2',mean_I2,...
  'mean_IN',mean_IN,...
  'mean_I2_over_H0',mean_I2_over_H0,...
  'mean_IN_over_H0',mean_IN_over_H0,...
  'mean_I2_over_IN',mean_I2_over_IN,...
  'mean_H1_over_H0',mean_H1_over_H0,...
  'mean_H2_over_H0',mean_H2_over_H0,...
  'std_H1',std_H1,...
  'std_H2',std_H2,...
  'std_H0',std_H0,...
  'std_I2',std_I2,...
  'std_IN',std_IN,...
  'std_I2_over_H0',std_I2_over_H0,...
  'std_IN_over_H0',std_IN_over_H0,...
  'std_I2_over_IN',std_I2_over_IN,...
  'std_H1_over_H0',std_H1_over_H0,...
  'std_H2_over_H0',std_H2_over_H0,...
  'perc_error_H1',perc_error_H1,...
  'perc_error_H2',perc_error_H2,...
  'perc_error_I2',perc_error_I2,...
  'perc_error_IN',perc_error_IN,...
  'perc_error_H0',perc_error_H0,...
  'perc_error_I2_over_H0',perc_error_I2_over_H0,...
  'perc_error_IN_over_H0',perc_error_IN_over_H0,...
  'perc_error_I2_over_IN',perc_error_I2_over_IN,...
  'perc_error_H1_over_H0',perc_error_H1_over_H0,...
  'perc_error_H2_over_H0',perc_error_H2_over_H0,...
  'report',report);

end

