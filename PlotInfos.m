function [] = PlotInfos(dataset,Ningreds,saveflag,titleflag,culture)
% Plots statistics related to entropies and informations over Ningreds for
% dataset
%
% Written by: DJ Strouse
% Last updated: Nov 18, 2013 by DJ Strouse
%
% INPUTS
% dataset [=] string in {'standard','meats','spices','cultures'}
% Ningreds [=] vector of positive integers = numbers of ingredients to use
% 
% (optional)
% saveflag [=] binary = indicates whether or not to save figure
% titleflag [=] binary = indicates whether or not to title the figures
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
% 
% OUTPUTS
% saves a figure
%
% FUTURE UPDATES
% 

% init
if ~exist('Ningreds','var')
  Ningreds = 4:2:18;
end
if ~exist('saveflag','var')
  saveflag = false; % don't save figs by default
end
if ~exist('titleflag','var')
  titleflag = true; % title figs by default
end
if strcmp(dataset,'cultures')
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
end

% load and process the data
H0s = zeros(length(Ningreds),1);
H1s = zeros(length(Ningreds),1);
H2s = zeros(length(Ningreds),1);
I2s = zeros(length(Ningreds),1);
INs = zeros(length(Ningreds),1);
I2oH0s = zeros(length(Ningreds),1);
INoH0s = zeros(length(Ningreds),1);
I2oINs = zeros(length(Ningreds),1);
H0_ebs = zeros(length(Ningreds),1);
H1_ebs = zeros(length(Ningreds),1);
H2_ebs = zeros(length(Ningreds),1);
I2_ebs = zeros(length(Ningreds),1);
IN_ebs = zeros(length(Ningreds),1);
I2oH0_ebs = zeros(length(Ningreds),1);
INoH0_ebs = zeros(length(Ningreds),1);
I2oIN_ebs = zeros(length(Ningreds),1);

for n = 1:length(Ningreds)
  N = Ningreds(n);
  if strcmp(dataset,'cultures')
    load(['data/cultures/',culture_str,'/',num2str(N),'ingred/large.mat']);
  else
    load(['data/',dataset,'/',num2str(N),'ingred/large.mat']);
  end
  H0s(n) = maxent1v2.mean_H0;
  H1s(n) = maxent1v2.mean_H1;
  H2s(n) = maxent1v2.mean_H2;
  I2s(n) = maxent1v2.mean_I2;
  INs(n) = maxent1v2.mean_IN;
  I2oH0s(n) = maxent1v2.mean_I2_over_H0;
  INoH0s(n) = maxent1v2.mean_IN_over_H0;
  I2oINs(n) = maxent1v2.mean_I2_over_IN;
  H0_ebs(n) = maxent1v2.std_H0;
  H1_ebs(n) = maxent1v2.std_H1;
  H2_ebs(n) = maxent1v2.std_H2;
  I2_ebs(n) = maxent1v2.std_I2;
  IN_ebs(n) = maxent1v2.std_IN;
  I2oH0_ebs(n) = maxent1v2.std_I2_over_H0;
  INoH0_ebs(n) = maxent1v2.std_IN_over_H0;
  I2oIN_ebs(n) = maxent1v2.std_I2_over_IN;
end
clear n;

% mkdir if necessary
if saveflag
  if strcmp(dataset,'cultures')
    if ~exist(['figures/cultures/',culture_str],'dir')
      mkdir(['figures/cultures/',culture_str]);
    end
  else
    if ~exist(['figures/',dataset],'dir')
      mkdir(['figures/',dataset]);
    end
  end
end

% make title string if necessary
if titleflag
  if strcmp(dataset,'cultures')
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
    title_str = [title_str,' (',num2str(Nrec_used),' recipes)'];
  elseif strcmp(dataset,'meats')
    title_str = 'Meats';
    title_str = [title_str,' (',num2str(Ningred),' total ingredients)'];
  elseif strcmp(dataset,'spices')
    title_str = 'Spices';
    title_str = [title_str,' (',num2str(Ningred),' total ingredients)'];
  elseif strcmp(dataset,'standard')
    title_str = 'Standard';
    % standard compared to both ingredient and recipe subsets so label with
    % both numbers
    title_str1 = [title_str,' (',num2str(Ningred),' total ingredients)'];
    title_str2 = [title_str,' (',num2str(Nrec),' recipes)'];
    clear title_str;
  end
end

% plot entropies
xfig = 10;
yfig = 10;
wfig = 30;
hfig = 20;
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(4));
shadedErrorBar(Ningreds,H1s,H1_ebs,{'Color',c(1,:),'LineWidth',2});
shadedErrorBar(Ningreds,H2s,H2_ebs,{'Color',c(2,:),'LineWidth',2});
shadedErrorBar(Ningreds,H0s,H0_ebs,{'Color',c(3,:),'LineWidth',2});
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
% ymax = max([H0s+H0_ebs;H1s+H1_ebs;H2s+H2_ebs]);
% ymax = 1.1*ymax;
ymax = 14;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('entropy')
if titleflag
  if strcmp(dataset,'standard')
    title(title_str1)
  else
    title(title_str)
  end
end
text(xmin+1,.95*ymax,'H1','FontSize',24,'Color',c(1,:))
text(xmin+1,.87*ymax,'H2','FontSize',24,'Color',c(2,:))
text(xmin+1,.79*ymax,'H0','FontSize',24,'Color',c(3,:))
if saveflag
  if strcmp(dataset,'cultures')
    file_name =...
      ['figures/cultures/',culture_str,'/ents_',num2str(max(Ningreds)),'in'];
  else
    file_name =...
      ['figures/',dataset,'/ents_',num2str(max(Ningreds)),'in'];
  end
  if titleflag
    file_name = [file_name,'_titled'];
    if strcmp(dataset,'standard')
      % need to further specify for standard since will also make a version
      % labeled with Nrec
      file_name = [file_name,'_wNingred'];
    end
  end
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
  if strcmp(dataset,'standard')&&titleflag
    title(title_str2)
    file_name =...
      ['figures/standard/ents_',num2str(max(Ningreds)),'in_titled_wNrec'];
    saveas(h,[file_name,'.fig']);
    export_fig(file_name,'-pdf','-transparent');
  end
end

% plot infos
xfig = 10;
yfig = 10;
wfig = 30;
hfig = 20;
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(4));
shadedErrorBar(Ningreds,I2s,I2_ebs,{'Color',c(1,:),'LineWidth',2});
shadedErrorBar(Ningreds,INs,IN_ebs,{'Color',c(2,:),'LineWidth',2});
shadedErrorBar(Ningreds,H0s,H0_ebs,{'Color',c(3,:),'LineWidth',2});
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
% ymax = max([H0s+H0_ebs;I2s+I2_ebs;INs+IN_ebs]);
% ymax = 1.1*ymax;
ymax = 11;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('entropy/information')
if titleflag
  if strcmp(dataset,'standard')
    title(title_str1)
  else
    title(title_str)
  end
end
text(xmin+1,.95*ymax,'H0','FontSize',24,'Color',c(3,:))
text(xmin+1,.87*ymax,'IN','FontSize',24,'Color',c(2,:))
text(xmin+1,.79*ymax,'I2','FontSize',24,'Color',c(1,:))
if saveflag
  if strcmp(dataset,'cultures')
    file_name =...
      ['figures/cultures/',culture_str,'/infos_',num2str(max(Ningreds)),'in'];
  else
    file_name =...
      ['figures/',dataset,'/infos_',num2str(max(Ningreds)),'in'];
  end
  if titleflag
    file_name = [file_name,'_titled'];
    if strcmp(dataset,'standard')
      % need to further specify for standard since will also make a version
      % labeled with Nrec
      file_name = [file_name,'_wNingred'];
    end
  end
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
  if strcmp(dataset,'standard')&&titleflag
    title(title_str2)
    file_name =...
      ['figures/standard/infos_',num2str(max(Ningreds)),'in_titled_wNrec'];
    saveas(h,[file_name,'.fig']);
    export_fig(file_name,'-pdf','-transparent');
  end
end

% plot info fractions
xfig = 10;
yfig = 10;
wfig = 30;
hfig = 20;
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(4));
shadedErrorBar(Ningreds,I2oH0s,I2oH0_ebs,{'Color',c(1,:),'LineWidth',2});
shadedErrorBar(Ningreds,INoH0s,INoH0_ebs,{'Color',c(2,:),'LineWidth',2});
shadedErrorBar(Ningreds,I2oINs,I2oIN_ebs,{'Color',c(3,:),'LineWidth',2});
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
% ymax = min(max([I2oH0s+I2oH0_ebs;INoH0s+INoH0_ebs;I2oINs+I2oIN_ebs]),1);
% ymax = 1.1*ymax;
ymax = 1;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('info fraction')
if titleflag
  if strcmp(dataset,'standard')
    title(title_str1)
  else
    title(title_str)
  end
end
text(xmin+1,.55*ymax,'I2/IN','FontSize',24,'Color',c(3,:))
text(xmin+1,.47*ymax,'IN/H0','FontSize',24,'Color',c(2,:))
text(xmin+1,.39*ymax,'I2/H0','FontSize',24,'Color',c(1,:))
if saveflag
  if strcmp(dataset,'cultures')
    file_name =...
      ['figures/cultures/',culture_str,'/infofracs_',num2str(max(Ningreds)),'in'];
  else
    file_name =...
      ['figures/',dataset,'/infofracs_',num2str(max(Ningreds)),'in'];
  end
  if titleflag
    file_name = [file_name,'_titled'];
    if strcmp(dataset,'standard')
      % need to further specify for standard since will also make a version
      % labeled with Nrec
      file_name = [file_name,'_wNingred'];
    end
  end
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
  if strcmp(dataset,'standard')&&titleflag
    title(title_str2)
    file_name =...
      ['figures/standard/infofracs_',num2str(max(Ningreds)),'in_titled_wNrec'];
    saveas(h,[file_name,'.fig']);
    export_fig(file_name,'-pdf','-transparent');
  end
end

end

