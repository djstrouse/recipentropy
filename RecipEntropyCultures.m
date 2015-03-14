function [] = RecipEntropyCultures(culture,Ningred_used,Nsubsamp)
% Fits a maxent model to recipes from a specified culture, using the
% Ningred_used most common ingredients and obtaining error bars by
% subsampling the data Nsubsamp times
%
% Written by: DJ Strouse
% Last updated: Aug 9, 2013 by DJ Strouse
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
% Ningred_used [=] positive integer = number of ingredients to use
% Nsubsamp [=] positive integer = number of data subsamplings (for ebs)
%
% OUTPUTS
% none

%% init 
addpath(genpath('/scratch/network/astrandb/recipentropy')) % TIGRESS
addpath(genpath('/scratch/network/dstrouse/recipentropy')) % TIGRESS
addpath(genpath('/scratch/network/kalindbl/recipentropy')) % TIGRESS
if ~exist('Nsubsamp')
  Nsubsamp = 25; % default number of subsamplings
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
diary(['recipentropy',culture_str,int2str(Ningred_used),'_diary.txt']);
load('data/recipes.mat');
Ningred = size(recipes_binary,2);
Nrec = size(recipes_binary,1);
if Ningred_used>Ningred
  error('Ningred_used>Ningred!')
end
plot_figs = false; % indicates whether to plot figures
save_figs = false; % indicates whether to save figures
save_data = true; % indicates whether to save data

%% use recipes of only the specified culture

culturemask = zeros(Nrec,1); % indicator function of culture identity
for r = 1:Nrec
  culturemask(r) = strcmp(recipes(r).region{1},culture_str);
end
clear r;
culturemask = culturemask==1; % convert to logical
recipes_cultural = recipes(culturemask);
recipes_binary_cultural = recipes_binary(culturemask,:);
Nrec_used = size(recipes_binary_cultural,1);
subsamp_size = round(.8*Nrec_used);

%% use Ningred_used most common ingredients only

[sorted_ingred_freq,i] = sort(ingred_freq,'descend');
to_use = i(1:Ningred_used); clear i;
recipes_binary_subsamp = recipes_binary_cultural(:,to_use)';
ingreds_subsamp = ingreds(to_use); clear to_use;
disp(sprintf('Using the following %i ingredients:',Ningred_used));
disp(ingreds_subsamp');

%% zipf plots (all)

if plot_figs

  xfig = 10;
  yfig = 10;
  wfig = 30;
  hfig = 20;
  buffer = .05;

  h = figure(1);
  set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
  plot(1:Ningred,sorted_ingred_freq,'.','LineWidth',2);
  set(gca,'XTick',1:Ningred_used)
  prettyplot
  xlim([0 Ningred_used+1])
  ylim([...
    sorted_ingred_freq(Ningred_used)-buffer*sorted_ingred_freq(Ningred_used)...
    sorted_ingred_freq(1)+buffer*sorted_ingred_freq(1)])
  xlabel('ingredient frequency rank')
  ylabel('ingredient frequency')
  if save_figs
    file_name = 'figures/zipf_linear_all';
    saveas(h,[file_name,'.fig']);
    export_fig(file_name,'-pdf','-transparent');
  end

  h = figure(2);
  set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
  plot(1:Ningred,sorted_ingred_freq./Nrec,'.','LineWidth',2);
  set(gca,'XTick',1:Ningred_used)
  prettyplot
  xlim([0 Ningred_used+1])
  ylim([...
    sorted_ingred_freq(Ningred_used)-buffer*sorted_ingred_freq(Ningred_used)...
    sorted_ingred_freq(1)+buffer*sorted_ingred_freq(1)]./Nrec)
  xlabel('ingredient frequency rank')
  ylabel('normalized ingredient frequency')
  if save_figs
    file_name = 'figures/zipf_linear_all_normalized';
    saveas(h,[file_name,'.fig']);
    export_fig(file_name,'-pdf','-transparent');
  end

  h = figure(3);
  set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
  semilogy(1:Ningred,sorted_ingred_freq,'.','LineWidth',2);
  set(gca,'XTick',1:Ningred_used)
  prettyplot
  xlim([0 Ningred_used+1])
  ylim([...
    sorted_ingred_freq(Ningred_used)-buffer*sorted_ingred_freq(Ningred_used)...
    sorted_ingred_freq(1)+buffer*sorted_ingred_freq(1)])
  xlabel('ingredient frequency rank')
  ylabel('ingredient frequency')
  if save_figs
    file_name = 'figures/zipf_semilog_all';
    saveas(h,[file_name,'.fig']);
    export_fig(file_name,'-pdf','-transparent');
  end

  h = figure(4);
  set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
  semilogy(1:Ningred,sorted_ingred_freq./Nrec,'.','LineWidth',2);
  set(gca,'XTick',1:Ningred_used)
  prettyplot
  xlim([0 Ningred_used+1])
  ylim([...
    sorted_ingred_freq(Ningred_used)-buffer*sorted_ingred_freq(Ningred_used)...
    sorted_ingred_freq(1)+buffer*sorted_ingred_freq(1)]./Nrec)
  xlabel('ingredient frequency rank')
  ylabel('normalized ingredient frequency')
  if save_figs
    file_name = 'figures/zipf_semilog_all_normalized';
    saveas(h,[file_name,'.fig']);
    export_fig(file_name,'-pdf','-transparent');
  end

  h = figure(5);
  set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
  loglog(1:Ningred,sorted_ingred_freq,'.','LineWidth',2);
  set(gca,'XTick',1:Ningred_used)
  prettyplot
  xlim([.9 Ningred_used+1])
  ylim([...
    sorted_ingred_freq(Ningred_used)-buffer*sorted_ingred_freq(Ningred_used)...
    sorted_ingred_freq(1)+buffer*sorted_ingred_freq(1)])
  xlabel('ingredient frequency rank')
  ylabel('ingredient frequency')
  if save_figs
    file_name = 'figures/zipf_loglog_all';
    saveas(h,[file_name,'.fig']);
    export_fig(file_name,'-pdf','-transparent');
  end

  h = figure(6);
  set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
  loglog(1:Ningred,sorted_ingred_freq./Nrec,'.','LineWidth',2);
  set(gca,'XTick',1:Ningred_used)
  prettyplot
  xlim([.9 Ningred_used+1])
  ylim([...
    sorted_ingred_freq(Ningred_used)-buffer*sorted_ingred_freq(Ningred_used)...
    sorted_ingred_freq(1)+buffer*sorted_ingred_freq(1)]./Nrec)
  xlabel('ingredient frequency rank')
  ylabel('normalized ingredient frequency')
  if save_figs
    file_name = 'figures/zipf_loglog_all_normalized';
    saveas(h,[file_name,'.fig']);
    export_fig(file_name,'-pdf','-transparent');
  end

end

%% labeled zipf plot (Ningred_used most common)

if plot_figs

  xfig = 10;
  yfig = 10;
  wfig = 30;
  hfig = 20;
  buffer = .05;

  h = figure(7);
  set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
  plot(1:Ningred_used,sorted_ingred_freq(1:Ningred_used),'.','LineWidth',2);
  set(gca,'XTick',1:Ningred_used,'XTickLabel',ingreds_subsamp)
  % ROTATE LABELS
  prettyplot
  xlim([0 Ningred_used+1])
  ylim([...
    sorted_ingred_freq(Ningred_used)-buffer*sorted_ingred_freq(Ningred_used)...
    sorted_ingred_freq(1)+buffer*sorted_ingred_freq(1)])
  xlabel('ingredient frequency rank')
  ylabel('ingredient frequency')
  if save_figs
    file_name = sprintf('figures/zipf_semilog_%iin',Ningred_used);
    saveas(h,[file_name,'.fig']);
    export_fig(file_name,'-pdf','-transparent');
  end

  h = figure(8);
  set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
  plot(1:Ningred_used,sorted_ingred_freq(1:Ningred_used)./Nrec,'.','LineWidth',2);
  set(gca,'XTick',1:Ningred_used,'XTickLabel',ingreds_subsamp)
  % ROTATE LABELS
  prettyplot
  xlim([0 Ningred_used+1])
  ylim([...
    sorted_ingred_freq(Ningred_used)-buffer*sorted_ingred_freq(Ningred_used)...
    sorted_ingred_freq(1)+buffer*sorted_ingred_freq(1)]./Nrec)
  xlabel('ingredient frequency rank')
  ylabel('normalized ingredient frequency')
  if save_figs
    file_name = sprintf('figures/zipf_semilog_%iin_normalized',Ningred_used);
    saveas(h,[file_name,'.fig']);
    export_fig(file_name,'-pdf','-transparent');
  end

  h = figure(9);
  set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
  semilogy(1:Ningred_used,sorted_ingred_freq(1:Ningred_used),'.','LineWidth',2);
  set(gca,'XTick',1:Ningred_used,'XTickLabel',ingreds_subsamp)
  % ROTATE LABELS
  prettyplot
  xlim([0 Ningred_used+1])
  ylim([...
    sorted_ingred_freq(Ningred_used)-buffer*sorted_ingred_freq(Ningred_used)...
    sorted_ingred_freq(1)+buffer*sorted_ingred_freq(1)])
  xlabel('ingredient frequency rank')
  ylabel('ingredient frequency')
  if save_figs
    file_name = sprintf('figures/zipf_semilog_%iin',Ningred_used);
    saveas(h,[file_name,'.fig']);
    export_fig(file_name,'-pdf','-transparent');
  end

  h = figure(10);
  set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
  semilogy(1:Ningred_used,sorted_ingred_freq(1:Ningred_used)./Nrec,'.','LineWidth',2);
  set(gca,'XTick',1:Ningred_used,'XTickLabel',ingreds_subsamp)
  % ROTATE LABELS
  prettyplot
  xlim([0 Ningred_used+1])
  ylim([...
    sorted_ingred_freq(Ningred_used)-buffer*sorted_ingred_freq(Ningred_used)...
    sorted_ingred_freq(1)+buffer*sorted_ingred_freq(1)]./Nrec)
  xlabel('ingredient frequency rank')
  ylabel('normalized ingredient frequency')
  if save_figs
    file_name = sprintf('figures/zipf_semilog_%iin_normalized',Ningred_used);
    saveas(h,[file_name,'.fig']);
    export_fig(file_name,'-pdf','-transparent');
  end

end

%% extract raw frequencies

disp('beginning subsampling raw freqs')
[ingred_freq_subsamp_small,ingred_freq_subsamp_large] =...
  FindFreqs(recipes_binary_subsamp,subsamp_size,Nsubsamp);
disp('finished subsampling raw freqs')

%% fit first order model

[maxent1_small,maxent1_large] =...
  FitMaxEnt(recipes_binary_subsamp,1,subsamp_size,Nsubsamp);

%% fit second order model

[maxent2_small,maxent2_large] =...
  FitMaxEnt(recipes_binary_subsamp,2,subsamp_size,Nsubsamp);

%% compare first and second order models

maxent1v2 = CompareMaxEnt12(...
  maxent1_large,maxent2_large,ingred_freq_subsamp_large,...
  'first-order maxent','second-order maxent','raw freqs');

%% save data
if save_data
  save(...
    ['data/recipentropy',culture_str,'_',int2str(Ningred_used),'in_small.mat'],...
    'maxent1_small',...
    'maxent2_small',...
    'ingred_freq_subsamp_small',...
    'maxent1v2',...
    'Ningred',...
    'Nrec',...
    'Nrec_used',...
    'Nsubsamp',...
    'subsamp_size',...
    'Ningred_used',...
    'culture',...
    'culture_str',...
    'recipes_cultural',...
    'recipes_binary_cultural',...
    '-v7.3');
  save(...
    ['data/recipentropy',culture_str,'_',int2str(Ningred_used),'in_large.mat'],...
    'maxent1_large',...
    'maxent2_large',...
    'ingred_freq_subsamp_large',...
    'maxent1v2',...
    'Ningred',...
    'Nrec',...
    'Nrec_used',...
    'Nsubsamp',...
    'subsamp_size',...
    'Ningred_used',...
    'ingreds',...
    'sorted_ingred_freq',...
    'ingreds_subsamp',...
    'recipes_binary',...
    'recipes_binary_subsamp',...
    'culture',...
    'culture_str',...
    'recipes_cultural',...
    'recipes_binary_cultural',...
    '-v7.3');
end
diary off;

%% fixes a problem when running this in cambridge
! rm -f nohup.out
  
end