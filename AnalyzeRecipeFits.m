%% init

close all; clear all;
save_figs = false; 
Ningreds = [4,6,8,10,12,14,16,20]';
Nfits = length(Ningreds);
maxent1v2s = cell(Nfits,1);
maxent1s = cell(Nfits,1);
maxent2s = cell(Nfits,1);
freqs = cell(Nfits,1);
load('data/recipentropy_4in_large.mat');
maxent1v2s{1} = maxent1v2;
maxent1s{1} = maxent1_large;
maxent2s{1} = maxent2_large;
freqs{1} = ingred_freq_subsamp_large;
clear ingred_freq_subsamp_large maxent1_large maxent2_large maxent1v2;
load('data/recipentropy_6in_large.mat');
maxent1v2s{2} = maxent1v2;
maxent1s{2} = maxent1_large;
maxent2s{2} = maxent2_large;
freqs{2} = ingred_freq_subsamp_large;
clear ingred_freq_subsamp_large maxent1_large maxent2_large maxent1v2;
load('data/recipentropy_8in_large.mat');
maxent1v2s{3} = maxent1v2;
maxent1s{3} = maxent1_large;
maxent2s{3} = maxent2_large;
freqs{3} = ingred_freq_subsamp_large;
clear ingred_freq_subsamp_large maxent1_large maxent2_large maxent1v2;
load('data/recipentropy_10in_large.mat');
maxent1v2s{4} = maxent1v2;
maxent1s{4} = maxent1_large;
maxent2s{4} = maxent2_large;
freqs{4} = ingred_freq_subsamp_large;
clear ingred_freq_subsamp_large maxent1_large maxent2_large maxent1v2;
load('data/recipentropy_12in_large.mat');
maxent1v2s{5} = maxent1v2;
maxent1s{5} = maxent1_large;
maxent2s{5} = maxent2_large;
freqs{5} = ingred_freq_subsamp_large;
clear ingred_freq_subsamp_large maxent1_large maxent2_large maxent1v2;
load('data/recipentropy_14in_large.mat');
maxent1v2s{6} = maxent1v2;
maxent1s{6} = maxent1_large;
maxent2s{6} = maxent2_large;
freqs{6} = ingred_freq_subsamp_large;
clear ingred_freq_subsamp_large maxent1_large maxent2_large maxent1v2;
load('data/recipentropy_16in_large.mat');
maxent1v2s{7} = maxent1v2;
maxent1s{7} = maxent1_large;
maxent2s{7} = maxent2_large;
freqs{7} = ingred_freq_subsamp_large;
clear ingred_freq_subsamp_large maxent1_large maxent2_large maxent1v2;
% load('data/recipentropy_18in_large.mat');
% maxent1v2s{8} = maxent1v2;
% maxent1s{8} = maxent1_large;
% maxent2s{8} = maxent2_large;
% freqs{8} = ingred_freq_subsamp_large;
% clear ingred_freq_subsamp_large maxent1_large maxent2_large maxent1v2;
load('data/recipentropy_20in_large.mat');
maxent1v2s{8} = maxent1v2;
maxent1s{8} = maxent1_large;
maxent2s{8} = maxent2_large;
freqs{8} = ingred_freq_subsamp_large;
clear ingred_freq_subsamp_large maxent1_large maxent2_large maxent1v2;

%% plot entropies

% put in vectors
H0s = zeros(Nfits,1);
H1s = zeros(Nfits,1);
H2s = zeros(Nfits,1);
H0_ebs = zeros(Nfits,1);
H1_ebs = zeros(Nfits,1);
H2_ebs = zeros(Nfits,1);
for i = 1:Nfits
  H0s(i) = maxent1v2s{i}.mean_H0;
  H1s(i) = maxent1v2s{i}.mean_H1;
  H2s(i) = maxent1v2s{i}.mean_H2;
  H0_ebs(i) = maxent1v2s{i}.std_H0;
  H1_ebs(i) = maxent1v2s{i}.std_H1;
  H2_ebs(i) = maxent1v2s{i}.std_H2;
end
clear i;

% plot it
h = figure(1);
hold off
hold on
c = colormap(pmkmp(4));
shadedErrorBar(Ningreds,H1s,H1_ebs,{'Color',c(1,:),'LineWidth',2});
shadedErrorBar(Ningreds,H2s,H2_ebs,{'Color',c(2,:),'LineWidth',2});
shadedErrorBar(Ningreds,H0s,H0_ebs,{'Color',c(3,:),'LineWidth',2});
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max([H0s+H0_ebs;H1s+H1_ebs;H2s+H2_ebs]);
ymax = 1.1*ymax;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('entropy')
text(xmin+1,.95*ymax,'H1','FontSize',24,'Color',c(1,:))
text(xmin+1,.87*ymax,'H2','FontSize',24,'Color',c(2,:))
text(xmin+1,.79*ymax,'H0','FontSize',24,'Color',c(3,:))
if save_figs
  file_name = sprintf('figures/entropies_%ito%iingreds',min(Ningreds),max(Ningreds));
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-eps','-transparent');
end

%% plot connected info and multi-info

% put in vectors
I2s = zeros(Nfits,1);
INs = zeros(Nfits,1);
I2_ebs = zeros(Nfits,1);
IN_ebs = zeros(Nfits,1);
for i = 1:Nfits
  I2s(i) = maxent1v2s{i}.mean_I2;
  INs(i) = maxent1v2s{i}.mean_IN;
  I2_ebs(i) = maxent1v2s{i}.std_I2;
  IN_ebs(i) = maxent1v2s{i}.std_IN;
end
clear i;

% plot it
h = figure(2);
hold off
hold on
c = colormap(pmkmp(4));
shadedErrorBar(Ningreds,H0s,H0_ebs,{'Color',c(1,:),'LineWidth',2});
shadedErrorBar(Ningreds,INs,IN_ebs,{'Color',c(2,:),'LineWidth',2});
shadedErrorBar(Ningreds,I2s,I2_ebs,{'Color',c(3,:),'LineWidth',2});
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max([H0s+H0_ebs;I2s+I2_ebs;INs+IN_ebs]);
ymax = 1.1*ymax;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('entropy/information')
text(xmin+1,.95*ymax,'H0','FontSize',24,'Color',c(1,:))
text(xmin+1,.87*ymax,'IN','FontSize',24,'Color',c(2,:))
text(xmin+1,.79*ymax,'I2','FontSize',24,'Color',c(3,:))
if save_figs
  file_name = sprintf('figures/informations_%ito%iingreds',min(Ningreds),max(Ningreds));
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-eps','-transparent');
end

%% plot info ratios

% put in vectors
I2_over_INs = zeros(Nfits,1);
IN_over_H0s = zeros(Nfits,1);
I2_over_H0s = zeros(Nfits,1);
H1_over_H0s = zeros(Nfits,1);
H2_over_H0s = zeros(Nfits,1);
I2_over_IN_ebs = zeros(Nfits,1);
IN_over_H0_ebs = zeros(Nfits,1);
I2_over_H0_ebs = zeros(Nfits,1);
H1_over_H0_ebs = zeros(Nfits,1);
H2_over_H0_ebs = zeros(Nfits,1);
for i = 1:Nfits
  I2_over_INs(i) = maxent1v2s{i}.mean_I2_over_IN;
  IN_over_H0s(i) = maxent1v2s{i}.mean_IN_over_H0;
  I2_over_H0s(i) = maxent1v2s{i}.mean_I2_over_H0;
  H1_over_H0s(i) = maxent1v2s{i}.mean_H1_over_H0;
  H2_over_H0s(i) = maxent1v2s{i}.mean_H2_over_H0;
  I2_over_IN_ebs(i) = maxent1v2s{i}.std_I2_over_IN;
  IN_over_H0_ebs(i) = maxent1v2s{i}.std_IN_over_H0;
  I2_over_H0_ebs(i) = maxent1v2s{i}.std_I2_over_H0;
  H1_over_H0_ebs(i) = maxent1v2s{i}.std_H1_over_H0;
  H2_over_H0_ebs(i) = maxent1v2s{i}.std_H2_over_H0;
end
clear i;

% plot it
h = figure(3);
hold off
hold on
c = colormap(pmkmp(6));
shadedErrorBar(Ningreds,I2_over_INs,I2_over_IN_ebs,{'Color',c(1,:),'LineWidth',2});
shadedErrorBar(Ningreds,IN_over_H0s,IN_over_H0_ebs,{'Color',c(2,:),'LineWidth',2});
shadedErrorBar(Ningreds,I2_over_H0s,I2_over_H0_ebs,{'Color',c(3,:),'LineWidth',2});
shadedErrorBar(Ningreds,H1_over_H0s,H1_over_H0_ebs,{'Color',c(4,:),'LineWidth',2});
shadedErrorBar(Ningreds,H2_over_H0s,H2_over_H0_ebs,{'Color',c(5,:),'LineWidth',2});
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = 1;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('fraction of entropy/information')
text(xmin+1,.58*ymax,'I2/IN','FontSize',24,'Color',c(1,:))
text(xmin+1,.5*ymax,'IN/H0','FontSize',24,'Color',c(2,:))
text(xmin+1,.42*ymax,'I2/H0','FontSize',24,'Color',c(3,:))
% add text!
%
if save_figs
  file_name = sprintf('figures/info_fractions_%ito%iingreds',min(Ningreds),max(Ningreds));
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-eps','-transparent');
end