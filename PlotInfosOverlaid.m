function [] = PlotInfosOverlaid(Ningreds,cultures,saveflag)
% Plots estatistics related to entropies and informations across Ningreds
% in separate figures each which include either (1) standard + cultures or
% (2) standard + meats + spices [4 figures total]
%
% Written by: DJ Strouse
% Last updated: Nov 22, 2013 by DJ Strouse
%
% INPUTS
% (optional)
% Ningreds [=] vector of positive integers = numbers of ingredients to use
% cultures [=] vector of positive integers = IDs of cultures included
%     (ranked by number of recs)
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
% saveflag [=] binary = indicates whether or not to save figure
% 
% OUTPUTS
% saves 16 figures
%
% FUTURE UPDATES
% 

% defaults
if ~exist('Ningreds','var')
  Ningreds = 4:2:18;
elseif isempty(Ningreds)
  Ningreds = 4:2:18;
end
if ~exist('cultures','var')
  cultures = 1:5;
elseif isempty(cultures)  
  cultures = 1:5;
end
if ~exist('saveflag','var')
  saveflag = false; % don't save figs by default
end
Ncultures = length(cultures);

% init
datasets = {'standard','meats','spices','cultures'};
Ntotal = length(datasets)+Ncultures-1;
allH0s = zeros(Ntotal,length(Ningreds));
allH1s = zeros(Ntotal,length(Ningreds));
allH2s = zeros(Ntotal,length(Ningreds));
allI2s = zeros(Ntotal,length(Ningreds));
allINs = zeros(Ntotal,length(Ningreds));
allI2oH0s = zeros(Ntotal,length(Ningreds));
allINoH0s = zeros(Ntotal,length(Ningreds));
allI2oINs = zeros(Ntotal,length(Ningreds));
allH0_ebs = zeros(Ntotal,length(Ningreds));
allH1_ebs = zeros(Ntotal,length(Ningreds));
allH2_ebs = zeros(Ntotal,length(Ningreds));
allI2_ebs = zeros(Ntotal,length(Ningreds));
allIN_ebs = zeros(Ntotal,length(Ningreds));
allI2oH0_ebs = zeros(Ntotal,length(Ningreds));
allINoH0_ebs = zeros(Ntotal,length(Ningreds));
allI2oIN_ebs = zeros(Ntotal,length(Ningreds));
allNingred_total = zeros(length(datasets)-1,1);
allNrec = zeros(Ncultures+1,1);

% iterate over datasets
for d = 1:length(datasets)
  dataset = datasets{d};
  % hacky way to iterate over all datasets
  if strcmp(dataset,'cultures')
    iters = cultures;
  else
    iters = 1;
  end
  for k = iters;
    if strcmp(dataset,'cultures')
      if k==1
        culture_str = 'NorthAmerican';
      elseif k==2
        culture_str = 'SouthernEuropean';
      elseif k==3
        culture_str = 'LatinAmerican';
      elseif k==4
        culture_str = 'WesternEuropean';
      elseif k==5
        culture_str = 'EastAsian';
      elseif k==6
        culture_str = 'MiddleEastern';
      elseif k==7
        culture_str = 'SouthAsian';
      elseif k==8
        culture_str = 'SoutheastAsian';
      elseif k==9
        culture_str = 'EasternEuropean';
      elseif k==10
        culture_str = 'African';
      elseif k==11
        culture_str = 'NorthernEuropean';
      end
      for n = 1:length(Ningreds)
        N = Ningreds(n);
        load(['data/cultures/',culture_str,'/',num2str(N),'ingred/large.mat']);
        allH0s(d+k-1,n) = maxent1v2.mean_H0;
        allH1s(d+k-1,n) = maxent1v2.mean_H1;
        allH2s(d+k-1,n) = maxent1v2.mean_H2;
        allI2s(d+k-1,n) = maxent1v2.mean_I2;
        allINs(d+k-1,n) = maxent1v2.mean_IN;
        allI2oH0s(d+k-1,n) = maxent1v2.mean_I2_over_H0;
        allINoH0s(d+k-1,n) = maxent1v2.mean_IN_over_H0;
        allI2oINs(d+k-1,n) = maxent1v2.mean_I2_over_IN;
        allH0_ebs(d+k-1,n) = maxent1v2.std_H0;
        allH1_ebs(d+k-1,n) = maxent1v2.std_H1;
        allH2_ebs(d+k-1,n) = maxent1v2.std_H2;
        allI2_ebs(d+k-1,n) = maxent1v2.std_I2;
        allIN_ebs(d+k-1,n) = maxent1v2.std_IN;
        allI2oH0_ebs(d+k-1,n) = maxent1v2.std_I2_over_H0;
        allINoH0_ebs(d+k-1,n) = maxent1v2.std_IN_over_H0;
        allI2oIN_ebs(d+k-1,n) = maxent1v2.std_I2_over_IN;
      end
      clear n;
      allNrec(k+1) = Nrec_used; clear Nrec_used;
    else % standard, meats, spices
      for n = 1:length(Ningreds)
        N = Ningreds(n);
        load(['data/',dataset,'/',num2str(N),'ingred/large.mat']);
        allH0s(d+k-1,n) = maxent1v2.mean_H0;
        allH1s(d+k-1,n) = maxent1v2.mean_H1;
        allH2s(d+k-1,n) = maxent1v2.mean_H2;
        allI2s(d+k-1,n) = maxent1v2.mean_I2;
        allINs(d+k-1,n) = maxent1v2.mean_IN;
        allI2oH0s(d+k-1,n) = maxent1v2.mean_I2_over_H0;
        allINoH0s(d+k-1,n) = maxent1v2.mean_IN_over_H0;
        allI2oINs(d+k-1,n) = maxent1v2.mean_I2_over_IN;
        allH0_ebs(d+k-1,n) = maxent1v2.std_H0;
        allH1_ebs(d+k-1,n) = maxent1v2.std_H1;
        allH2_ebs(d+k-1,n) = maxent1v2.std_H2;
        allI2_ebs(d+k-1,n) = maxent1v2.std_I2;
        allIN_ebs(d+k-1,n) = maxent1v2.std_IN;
        allI2oH0_ebs(d+k-1,n) = maxent1v2.std_I2_over_H0;
        allINoH0_ebs(d+k-1,n) = maxent1v2.std_IN_over_H0;
        allI2oIN_ebs(d+k-1,n) = maxent1v2.std_I2_over_IN;
      end
      clear n;
      if strcmp(dataset,'standard')
        allNrec(1) = Nrec; clear Nrec;
      end
      allNingred_total(d) = Ningred; clear Ningred;
    end
  end
  clear k;
end
clear d;

% mkdir if necessary
if saveflag
  if ~exist('figures/combined/','dir')
    mkdir('figures/combined/');
  end
end

% PLOTS
xfig = 10;
yfig = 10;
wfig = 30;
hfig = 20;

% plot H0 across ingredient subsets
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(4));
shadedErrorBar(Ningreds,allH0s(1,:)',allH0_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
shadedErrorBar(Ningreds,allH0s(2,:)',allH0_ebs(2,:)',{'Color',c(2,:),'LineWidth',2});
shadedErrorBar(Ningreds,allH0s(3,:)',allH0_ebs(3,:)',{'Color',c(3,:),'LineWidth',2});
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(...
  [allH0s(1,:)+allH0_ebs(1,:),...
  allH0s(2,:)+allH0_ebs(2,:),...
  allH0s(3,:)+allH0_ebs(3,:)]);
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('H_{0}')
text(...
  xmin+.5,.95*ymax,...
  ['standard (',num2str(allNingred_total(1)),' ingredients)'],...
  'FontSize',24,'Color',c(1,:))
text(...
  xmin+.5,.87*ymax,...
  ['spices (',num2str(allNingred_total(3)),' ingredients)'],...
  'FontSize',24,'Color',c(3,:))
text(...
  xmin+.5,.79*ymax,...
  ['meats (',num2str(allNingred_total(2)),' ingredients)'],...
  'FontSize',24,'Color',c(2,:))
if saveflag
  file_name = ['figures/combined/H0_insubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

% plot H1 across ingredient subsets
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(4));
shadedErrorBar(Ningreds,allH1s(1,:)',allH1_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
shadedErrorBar(Ningreds,allH1s(2,:)',allH1_ebs(2,:)',{'Color',c(2,:),'LineWidth',2});
shadedErrorBar(Ningreds,allH1s(3,:)',allH1_ebs(3,:)',{'Color',c(3,:),'LineWidth',2});
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(...
  [allH1s(1,:)+allH1_ebs(1,:),...
  allH1s(2,:)+allH1_ebs(2,:),...
  allH1s(3,:)+allH1_ebs(3,:)]);
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('H_{1}')
text(...
  xmin+.5,.95*ymax,...
  ['standard (',num2str(allNingred_total(1)),' ingredients)'],...
  'FontSize',24,'Color',c(1,:))
text(...
  xmin+.5,.87*ymax,...
  ['spices (',num2str(allNingred_total(3)),' ingredients)'],...
  'FontSize',24,'Color',c(3,:))
text(...
  xmin+.5,.79*ymax,...
  ['meats (',num2str(allNingred_total(2)),' ingredients)'],...
  'FontSize',24,'Color',c(2,:))
if saveflag
  file_name = ['figures/combined/H1_insubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

% plot H2 across ingredient subsets
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(4));
shadedErrorBar(Ningreds,allH2s(1,:)',allH2_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
shadedErrorBar(Ningreds,allH2s(2,:)',allH2_ebs(2,:)',{'Color',c(2,:),'LineWidth',2});
shadedErrorBar(Ningreds,allH2s(3,:)',allH2_ebs(3,:)',{'Color',c(3,:),'LineWidth',2});
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(...
  [allH2s(1,:)+allH2_ebs(1,:),...
  allH2s(2,:)+allH2_ebs(2,:),...
  allH2s(3,:)+allH2_ebs(3,:)]);
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('H_{2}')
text(...
  xmin+.5,.95*ymax,...
  ['standard (',num2str(allNingred_total(1)),' ingredients)'],...
  'FontSize',24,'Color',c(1,:))
text(...
  xmin+.5,.87*ymax,...
  ['spices (',num2str(allNingred_total(3)),' ingredients)'],...
  'FontSize',24,'Color',c(3,:))
text(...
  xmin+.5,.79*ymax,...
  ['meats (',num2str(allNingred_total(2)),' ingredients)'],...
  'FontSize',24,'Color',c(2,:))
if saveflag
  file_name = ['figures/combined/H2_insubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

% plot I2 across ingredient subsets
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(4));
shadedErrorBar(Ningreds,allI2s(1,:)',allI2_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
shadedErrorBar(Ningreds,allI2s(2,:)',allI2_ebs(2,:)',{'Color',c(2,:),'LineWidth',2});
shadedErrorBar(Ningreds,allI2s(3,:)',allI2_ebs(3,:)',{'Color',c(3,:),'LineWidth',2});
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(...
  [allI2s(1,:)+allI2_ebs(1,:),...
  allI2s(2,:)+allI2_ebs(2,:),...
  allI2s(3,:)+allI2_ebs(3,:)]);
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('I_{2}')
text(...
  xmin+.5,.95*ymax,...
  ['standard (',num2str(allNingred_total(1)),' ingredients)'],...
  'FontSize',24,'Color',c(1,:))
text(...
  xmin+.5,.87*ymax,...
  ['spices (',num2str(allNingred_total(3)),' ingredients)'],...
  'FontSize',24,'Color',c(3,:))
text(...
  xmin+.5,.79*ymax,...
  ['meats (',num2str(allNingred_total(2)),' ingredients)'],...
  'FontSize',24,'Color',c(2,:))
if saveflag
  file_name = ['figures/combined/I2_insubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

% plot IN across ingredient subsets
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(4));
shadedErrorBar(Ningreds,allINs(1,:)',allIN_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
shadedErrorBar(Ningreds,allINs(2,:)',allIN_ebs(2,:)',{'Color',c(2,:),'LineWidth',2});
shadedErrorBar(Ningreds,allINs(3,:)',allIN_ebs(3,:)',{'Color',c(3,:),'LineWidth',2});
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(...
  [allINs(1,:)+allIN_ebs(1,:),...
  allINs(2,:)+allIN_ebs(2,:),...
  allINs(3,:)+allIN_ebs(3,:)]);
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('I_{N}')
text(...
  xmin+.5,.95*ymax,...
  ['standard (',num2str(allNingred_total(1)),' ingredients)'],...
  'FontSize',24,'Color',c(1,:))
text(...
  xmin+.5,.87*ymax,...
  ['spices (',num2str(allNingred_total(3)),' ingredients)'],...
  'FontSize',24,'Color',c(3,:))
text(...
  xmin+.5,.79*ymax,...
  ['meats (',num2str(allNingred_total(2)),' ingredients)'],...
  'FontSize',24,'Color',c(2,:))
if saveflag
  file_name = ['figures/combined/IN_insubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

% plot I2/H0 across ingredient subsets
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(4));
shadedErrorBar(Ningreds,allI2oH0s(1,:)',allI2oH0_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
shadedErrorBar(Ningreds,allI2oH0s(2,:)',allI2oH0_ebs(2,:)',{'Color',c(2,:),'LineWidth',2});
shadedErrorBar(Ningreds,allI2oH0s(3,:)',allI2oH0_ebs(3,:)',{'Color',c(3,:),'LineWidth',2});
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(...
  [allI2oH0s(1,:)+allI2oH0_ebs(1,:),...
  allI2oH0s(2,:)+allI2oH0_ebs(2,:),...
  allI2oH0s(3,:)+allI2oH0_ebs(3,:)]);
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 min(ymax,1)])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('I_{2}/H_{0}')
text(...
  xmin+.5,.95*ymax,...
  ['standard (',num2str(allNingred_total(1)),' ingredients)'],...
  'FontSize',24,'Color',c(1,:))
text(...
  xmin+.5,.87*ymax,...
  ['spices (',num2str(allNingred_total(3)),' ingredients)'],...
  'FontSize',24,'Color',c(3,:))
text(...
  xmin+.5,.79*ymax,...
  ['meats (',num2str(allNingred_total(2)),' ingredients)'],...
  'FontSize',24,'Color',c(2,:))
if saveflag
  file_name = ['figures/combined/I2oH0_insubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

% plot IN/H0 across ingredient subsets
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(4));
shadedErrorBar(Ningreds,allINoH0s(1,:)',allINoH0_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
shadedErrorBar(Ningreds,allINoH0s(2,:)',allINoH0_ebs(2,:)',{'Color',c(2,:),'LineWidth',2});
shadedErrorBar(Ningreds,allINoH0s(3,:)',allINoH0_ebs(3,:)',{'Color',c(3,:),'LineWidth',2});
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(...
  [allINoH0s(1,:)+allINoH0_ebs(1,:),...
  allINoH0s(2,:)+allINoH0_ebs(2,:),...
  allINoH0s(3,:)+allINoH0_ebs(3,:)]);
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 min(ymax,1)])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('I_{N}/H_{0}')
text(...
  xmin+.5,.95*ymax,...
  ['standard (',num2str(allNingred_total(1)),' ingredients)'],...
  'FontSize',24,'Color',c(1,:))
text(...
  xmin+.5,.87*ymax,...
  ['spices (',num2str(allNingred_total(3)),' ingredients)'],...
  'FontSize',24,'Color',c(3,:))
text(...
  xmin+.5,.79*ymax,...
  ['meats (',num2str(allNingred_total(2)),' ingredients)'],...
  'FontSize',24,'Color',c(2,:))
if saveflag
  file_name = ['figures/combined/INoH0_insubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

% plot I2/IN across ingredient subsets
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(4));
shadedErrorBar(Ningreds,allI2oINs(1,:)',allI2oIN_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
shadedErrorBar(Ningreds,allI2oINs(2,:)',allI2oIN_ebs(2,:)',{'Color',c(2,:),'LineWidth',2});
shadedErrorBar(Ningreds,allI2oINs(3,:)',allI2oIN_ebs(3,:)',{'Color',c(3,:),'LineWidth',2});
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(...
  [allI2oINs(1,:)+allI2oIN_ebs(1,:),...
  allI2oINs(2,:)+allI2oIN_ebs(2,:),...
  allI2oINs(3,:)+allI2oIN_ebs(3,:)]);
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 min(ymax,1)])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('I_{2}/I_{N}')
text(...
  xmin+.5,.3*ymax,...
  ['standard (',num2str(allNingred_total(1)),' ingredients)'],...
  'FontSize',24,'Color',c(1,:))
text(...
  xmin+.5,.22*ymax,...
  ['spices (',num2str(allNingred_total(3)),' ingredients)'],...
  'FontSize',24,'Color',c(3,:))
text(...
  xmin+.5,.14*ymax,...
  ['meats (',num2str(allNingred_total(2)),' ingredients)'],...
  'FontSize',24,'Color',c(2,:))
if saveflag
  file_name = ['figures/combined/I2oIN_insubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

% plot H0 across recipe subsets
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(Ncultures+3));
shadedErrorBar(Ningreds,allH0s(1,:)',allH0_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
for i = 1:Ncultures
  shadedErrorBar(Ningreds,allH0s(3+i,:)',allH0_ebs(3+i,:)',{'Color',c(1+i,:),'LineWidth',2});
end
clear i;
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(allH0s(1,:)+allH0_ebs(1,:));
for i = 1:Ncultures
  ymax = max(ymax,max(allH0s(3+i,:)+allH0_ebs(3+i,:)));
end
clear i;
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('H_{0}')
text(...
  xmin+.5,.95*ymax,...
  ['standard (',num2str(allNingred_total(1)),' ingredients)'],...
  'FontSize',18,'Color',c(1,:))
for i = 1:Ncultures
  culture = cultures(i); % culture index
  if culture==1
    title_str = 'North American';
  elseif culture==2
    title_str = 'Southern European';
  elseif culture==3
    title_str = 'Latin American';
  elseif culture==4
    title_str = 'Western European';
  elseif culture==5
    title_str = 'East Asian';
  elseif culture==6
    title_str = 'Middle Eastern';
  elseif culture==7
    title_str = 'South Asian';
  elseif culture==8
    title_str = 'Southeast Asian';
  elseif culture==9
    title_str = 'Eastern European';
  elseif culture==10
    title_str = 'African';
  elseif culture==11
    title_str = 'Northern European';
  end
  text(...
    xmin+.5,(.95-.055*i)*ymax,...
    [title_str,' (',num2str(allNrec(1+i)),' recipes)'],...
    'FontSize',18,'Color',c(1+i,:))
end
clear i culture title_str;
if saveflag
  file_name = ['figures/combined/H0_recsubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

% plot H1 across recipe subsets
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(Ncultures+3));
shadedErrorBar(Ningreds,allH1s(1,:)',allH1_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
for i = 1:Ncultures
  shadedErrorBar(Ningreds,allH1s(3+i,:)',allH1_ebs(3+i,:)',{'Color',c(1+i,:),'LineWidth',2});
end
clear i;
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(allH1s(1,:)+allH1_ebs(1,:));
for i = 1:Ncultures
  ymax = max(ymax,max(allH1s(3+i,:)+allH1_ebs(3+i,:)));
end
clear i;
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('H_{1}')
text(...
  xmin+.5,.95*ymax,...
  ['standard (',num2str(allNingred_total(1)),' ingredients)'],...
  'FontSize',18,'Color',c(1,:))
for i = 1:Ncultures
  culture = cultures(i); % culture index
  if culture==1
    title_str = 'North American';
  elseif culture==2
    title_str = 'Southern European';
  elseif culture==3
    title_str = 'Latin American';
  elseif culture==4
    title_str = 'Western European';
  elseif culture==5
    title_str = 'East Asian';
  elseif culture==6
    title_str = 'Middle Eastern';
  elseif culture==7
    title_str = 'South Asian';
  elseif culture==8
    title_str = 'Southeast Asian';
  elseif culture==9
    title_str = 'Eastern European';
  elseif culture==10
    title_str = 'African';
  elseif culture==11
    title_str = 'Northern European';
  end
  text(...
    xmin+.5,(.95-.055*i)*ymax,...
    [title_str,' (',num2str(allNrec(1+i)),' recipes)'],...
    'FontSize',18,'Color',c(1+i,:))
end
clear i culture title_str;
if saveflag
  file_name = ['figures/combined/H1_recsubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

% plot H2 across recipe subsets
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(Ncultures+3));
shadedErrorBar(Ningreds,allH2s(1,:)',allH2_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
for i = 1:Ncultures
  shadedErrorBar(Ningreds,allH2s(3+i,:)',allH2_ebs(3+i,:)',{'Color',c(1+i,:),'LineWidth',2});
end
clear i;
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(allH2s(1,:)+allH2_ebs(1,:));
for i = 1:Ncultures
  ymax = max(ymax,max(allH2s(3+i,:)+allH2_ebs(3+i,:)));
end
clear i;
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('H_{2}')
text(...
  xmin+.5,.95*ymax,...
  ['standard (',num2str(allNingred_total(1)),' ingredients)'],...
  'FontSize',18,'Color',c(1,:))
for i = 1:Ncultures
  culture = cultures(i); % culture index
  if culture==1
    title_str = 'North American';
  elseif culture==2
    title_str = 'Southern European';
  elseif culture==3
    title_str = 'Latin American';
  elseif culture==4
    title_str = 'Western European';
  elseif culture==5
    title_str = 'East Asian';
  elseif culture==6
    title_str = 'Middle Eastern';
  elseif culture==7
    title_str = 'South Asian';
  elseif culture==8
    title_str = 'Southeast Asian';
  elseif culture==9
    title_str = 'Eastern European';
  elseif culture==10
    title_str = 'African';
  elseif culture==11
    title_str = 'Northern European';
  end
  text(...
    xmin+.5,(.95-.055*i)*ymax,...
    [title_str,' (',num2str(allNrec(1+i)),' recipes)'],...
    'FontSize',18,'Color',c(1+i,:))
end
clear i culture title_str;
if saveflag
  file_name = ['figures/combined/H2_recsubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

% plot I2 across recipe subsets
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(Ncultures+3));
shadedErrorBar(Ningreds,allI2s(1,:)',allI2_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
for i = 1:Ncultures
  shadedErrorBar(Ningreds,allI2s(3+i,:)',allI2_ebs(3+i,:)',{'Color',c(1+i,:),'LineWidth',2});
end
clear i;
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(allI2s(1,:)+allI2_ebs(1,:));
for i = 1:Ncultures
  ymax = max(ymax,max(allI2s(3+i,:)+allI2_ebs(3+i,:)));
end
clear i;
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('I_{2}')
text(...
  xmin+.5,.95*ymax,...
  ['standard (',num2str(allNingred_total(1)),' ingredients)'],...
  'FontSize',18,'Color',c(1,:))
for i = 1:Ncultures
  culture = cultures(i); % culture index
  if culture==1
    title_str = 'North American';
  elseif culture==2
    title_str = 'Southern European';
  elseif culture==3
    title_str = 'Latin American';
  elseif culture==4
    title_str = 'Western European';
  elseif culture==5
    title_str = 'East Asian';
  elseif culture==6
    title_str = 'Middle Eastern';
  elseif culture==7
    title_str = 'South Asian';
  elseif culture==8
    title_str = 'Southeast Asian';
  elseif culture==9
    title_str = 'Eastern European';
  elseif culture==10
    title_str = 'African';
  elseif culture==11
    title_str = 'Northern European';
  end
  text(...
    xmin+.5,(.95-.055*i)*ymax,...
    [title_str,' (',num2str(allNrec(1+i)),' recipes)'],...
    'FontSize',18,'Color',c(1+i,:))
end
clear i culture title_str;
if saveflag
  file_name = ['figures/combined/I2_recsubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

% plot IN across recipe subsets
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(Ncultures+3));
shadedErrorBar(Ningreds,allINs(1,:)',allIN_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
for i = 1:Ncultures
  shadedErrorBar(Ningreds,allINs(3+i,:)',allIN_ebs(3+i,:)',{'Color',c(1+i,:),'LineWidth',2});
end
clear i;
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(allINs(1,:)+allIN_ebs(1,:));
for i = 1:Ncultures
  ymax = max(ymax,max(allINs(3+i,:)+allIN_ebs(3+i,:)));
end
clear i;
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('I_{N}')
text(...
  xmin+.5,.95*ymax,...
  ['standard (',num2str(allNingred_total(1)),' ingredients)'],...
  'FontSize',18,'Color',c(1,:))
for i = 1:Ncultures
  culture = cultures(i); % culture index
  if culture==1
    title_str = 'North American';
  elseif culture==2
    title_str = 'Southern European';
  elseif culture==3
    title_str = 'Latin American';
  elseif culture==4
    title_str = 'Western European';
  elseif culture==5
    title_str = 'East Asian';
  elseif culture==6
    title_str = 'Middle Eastern';
  elseif culture==7
    title_str = 'South Asian';
  elseif culture==8
    title_str = 'Southeast Asian';
  elseif culture==9
    title_str = 'Eastern European';
  elseif culture==10
    title_str = 'African';
  elseif culture==11
    title_str = 'Northern European';
  end
  text(...
    xmin+.5,(.95-.055*i)*ymax,...
    [title_str,' (',num2str(allNrec(1+i)),' recipes)'],...
    'FontSize',18,'Color',c(1+i,:))
end
clear i culture title_str;
if saveflag
  file_name = ['figures/combined/IN_recsubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

% plot I2/H0 across recipe subsets
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(Ncultures+3));
shadedErrorBar(Ningreds,allI2oH0s(1,:)',allI2oH0_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
for i = 1:Ncultures
  shadedErrorBar(Ningreds,allI2oH0s(3+i,:)',allI2oH0_ebs(3+i,:)',{'Color',c(1+i,:),'LineWidth',2});
end
clear i;
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(allI2oH0s(1,:)+allI2oH0_ebs(1,:));
for i = 1:Ncultures
  ymax = max(ymax,max(allI2oH0s(3+i,:)+allI2oH0_ebs(3+i,:)));
end
clear i;
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 min(ymax,1)])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('I_{2}/H_{0}')
text(...
  xmin+.5,.95*ymax,...
  ['standard (',num2str(allNingred_total(1)),' ingredients)'],...
  'FontSize',18,'Color',c(1,:))
for i = 1:Ncultures
  culture = cultures(i); % culture index
  if culture==1
    title_str = 'North American';
  elseif culture==2
    title_str = 'Southern European';
  elseif culture==3
    title_str = 'Latin American';
  elseif culture==4
    title_str = 'Western European';
  elseif culture==5
    title_str = 'East Asian';
  elseif culture==6
    title_str = 'Middle Eastern';
  elseif culture==7
    title_str = 'South Asian';
  elseif culture==8
    title_str = 'Southeast Asian';
  elseif culture==9
    title_str = 'Eastern European';
  elseif culture==10
    title_str = 'African';
  elseif culture==11
    title_str = 'Northern European';
  end
  text(...
    xmin+.5,(.95-.055*i)*ymax,...
    [title_str,' (',num2str(allNrec(1+i)),' recipes)'],...
    'FontSize',18,'Color',c(1+i,:))
end
clear i culture title_str;
if saveflag
  file_name = ['figures/combined/I2oH0_recsubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

% plot IN/H0 across recipe subsets
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(Ncultures+3));
shadedErrorBar(Ningreds,allINoH0s(1,:)',allINoH0_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
for i = 1:Ncultures
  shadedErrorBar(Ningreds,allINoH0s(3+i,:)',allINoH0_ebs(3+i,:)',{'Color',c(1+i,:),'LineWidth',2});
end
clear i;
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(allINoH0s(1,:)+allINoH0_ebs(1,:));
for i = 1:Ncultures
  ymax = max(ymax,max(allINoH0s(3+i,:)+allINoH0_ebs(3+i,:)));
end
clear i;
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 min(ymax,1)])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('I_{N}/H_{0}')
text(...
  xmin+.5,.95*ymax,...
  ['standard (',num2str(allNingred_total(1)),' ingredients)'],...
  'FontSize',18,'Color',c(1,:))
for i = 1:Ncultures
  culture = cultures(i); % culture index
  if culture==1
    title_str = 'North American';
  elseif culture==2
    title_str = 'Southern European';
  elseif culture==3
    title_str = 'Latin American';
  elseif culture==4
    title_str = 'Western European';
  elseif culture==5
    title_str = 'East Asian';
  elseif culture==6
    title_str = 'Middle Eastern';
  elseif culture==7
    title_str = 'South Asian';
  elseif culture==8
    title_str = 'Southeast Asian';
  elseif culture==9
    title_str = 'Eastern European';
  elseif culture==10
    title_str = 'African';
  elseif culture==11
    title_str = 'Northern European';
  end
  text(...
    xmin+.5,(.95-.055*i)*ymax,...
    [title_str,' (',num2str(allNrec(1+i)),' recipes)'],...
    'FontSize',18,'Color',c(1+i,:))
end
clear i culture title_str;
if saveflag
  file_name = ['figures/combined/INoH0_recsubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

% plot I2/IN across recipe subsets
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(Ncultures+3));
shadedErrorBar(Ningreds,allI2oINs(1,:)',allI2oIN_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
for i = 1:Ncultures
  shadedErrorBar(Ningreds,allI2oINs(3+i,:)',allI2oIN_ebs(3+i,:)',{'Color',c(1+i,:),'LineWidth',2});
end
clear i;
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(allI2oINs(1,:)+allI2oIN_ebs(1,:));
for i = 1:Ncultures
  ymax = max(ymax,max(allI2oINs(3+i,:)+allI2oIN_ebs(3+i,:)));
end
clear i;
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 min(ymax,1)])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('I_{2}/I_{N}')
text(...
  xmin+.5,.37*ymax,...
  ['standard (',num2str(allNrec(1)),' recipes)'],...
  'FontSize',18,'Color',c(1,:))
for i = 1:Ncultures
  culture = cultures(i); % culture index
  if culture==1
    title_str = 'North American';
  elseif culture==2
    title_str = 'Southern European';
  elseif culture==3
    title_str = 'Latin American';
  elseif culture==4
    title_str = 'Western European';
  elseif culture==5
    title_str = 'East Asian';
  elseif culture==6
    title_str = 'Middle Eastern';
  elseif culture==7
    title_str = 'South Asian';
  elseif culture==8
    title_str = 'Southeast Asian';
  elseif culture==9
    title_str = 'Eastern European';
  elseif culture==10
    title_str = 'African';
  elseif culture==11
    title_str = 'Northern European';
  end
  text(...
    xmin+.5,(.37-.06*i)*ymax,...
    [title_str,' (',num2str(allNrec(1+i)),' recipes)'],...
    'FontSize',18,'Color',c(1+i,:))
end
clear i culture title_str;
if saveflag
  file_name = ['figures/combined/I2oIN_recsubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

end

