%% init

close all; clear all;
Ningred_used = 10;
Nsubsamp = 25;
diary(sprintf('recipentropyspices%i_diary.txt',Ningred_used));
load('data/spices.mat');
recipes_binary = final_spices_binary; clear final_spices_binary;
ingreds = final_spices; clear final_spices;
Ningred = size(recipes_binary,2);
Nrec = size(recipes_binary,1);
subsamp_size = round(.8*Nrec);
if Ningred_used>Ningred
  error('Ningred_used>Ningred!')
end
plot_figs = false; % indicates whether to plot figures
save_figs = false; % indicates whether to save figures
save_data = true; % indicates whether to save data

%% use Ningred_used most common ingredients only

ingred_freq = sum(recipes_binary,1);
[sorted_ingred_freq,i] = sort(ingred_freq,'descend');
to_use = i(1:Ningred_used); clear i;
recipes_binary_subsamp = recipes_binary(:,to_use)';
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
    sprintf('data/recipentropyspices_%iin_small.mat',Ningred_used),...
    'maxent1_small',...
    'maxent2_small',...
    'ingred_freq_subsamp_small',...
    'maxent1v2',...
    'Ningred',...
    'Nrec',...
    'Nsubsamp',...
    'subsamp_size',...
    'Ningred_used',...
    '-v7.3');
  save(...
    sprintf('data/recipentropyspices_%iin_large.mat',Ningred_used),...
    'maxent1_large',...
    'maxent2_large',...
    'ingred_freq_subsamp_large',...
    'maxent1v2',...
    'Ningred',...
    'Nrec',...
    'Nsubsamp',...
    'subsamp_size',...
    'Ningred_used',...
    'ingreds',...
    'sorted_ingred_freq',...
    'ingreds_subsamp',...
    'recipes_binary',...
    'recipes_binary_subsamp',...
    '-v7.3');
end
diary off;

%% fixes a problem when running this in cambridge
! rm -f nohup.out