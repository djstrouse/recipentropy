function [] = PlotJSOverlaid(Ningreds,cultures,saveflag)
% Plots JS divergences between P0/P1, P0/P2, and P1/P2 across Ningreds in
% separate figures each which include either (1) standard + cultures or (2)
% standard + meats + spices [6 figures total]
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
% saves 6 figures
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
allJS01s = zeros(Ntotal,length(Ningreds));
allJS02s = zeros(Ntotal,length(Ningreds));
allJS12s = zeros(Ntotal,length(Ningreds));
allJS01_ebs = zeros(Ntotal,length(Ningreds));
allJS02_ebs = zeros(Ntotal,length(Ningreds));
allJS12_ebs = zeros(Ntotal,length(Ningreds));
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
    % check if JS calculations have already been made
    docalc = true;
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
      % first check if dir exists
      if exist(['data/cultures/',culture_str,'/JS_',num2str(max(Ningreds)),'in.mat'],'file')
        % if it does, load it
        load(['data/cultures/',culture_str,'/JS_',num2str(max(Ningreds)),'in.mat']);
        % make sure same Ningreds were used previously
        if isequal(Ningreds,Ningreds_JS)
          % if so, no need to recalc JSs
          docalc = false;
          % save everything to the larger data structure
          allJS01s(d+k-1,:) = JS01s';
          allJS02s(d+k-1,:) = JS02s';
          allJS12s(d+k-1,:) = JS12s';
          allJS01_ebs(d+k-1,:) = JS01_ebs';
          allJS02_ebs(d+k-1,:) = JS02_ebs';
          allJS12_ebs(d+k-1,:) = JS12_ebs';
          allNrec(k+1) = Nrec;
          clear JS01s JS02s JS12s JS01_ebs JS02_ebs JS12_ebs;
          clear Ningred_total Ningreds_JS Nrec;
        end
      end
    else
      if exist(['data/',dataset,'/JS_',num2str(max(Ningreds)),'in.mat'],'file')
        load(['data/',dataset,'/JS_',num2str(max(Ningreds)),'in.mat']);
        if isequal(Ningreds,Ningreds_JS)
          docalc = false;
          allJS01s(d+k-1,:) = JS01s';
          allJS02s(d+k-1,:) = JS02s';
          allJS12s(d+k-1,:) = JS12s';
          allJS01_ebs(d+k-1,:) = JS01_ebs';
          allJS02_ebs(d+k-1,:) = JS02_ebs';
          allJS12_ebs(d+k-1,:) = JS12_ebs';
          allNingred_total(d) = Ningred_total;
          if strcmp(dataset,'standard')
            allNrec(1) = Nrec;
          end
          clear JS01s JS02s JS12s JS01_ebs JS02_ebs JS12_ebs;
          clear Ningred_total Ningreds_JS Nrec;
        end
      end
    end

    % load and process the data
    if docalc
      JS01s = zeros(length(Ningreds),1);
      JS01_ebs = zeros(length(Ningreds),1);
      JS02s = zeros(length(Ningreds),1);
      JS02_ebs = zeros(length(Ningreds),1);
      JS12s = zeros(length(Ningreds),1);
      JS12_ebs = zeros(length(Ningreds),1);
      for n = 1:length(Ningreds)
        N = Ningreds(n);
        if strcmp(dataset,'cultures')
          load(['data/cultures/',culture_str,'/',num2str(N),'ingred/large.mat']);
        else
          load(['data/',dataset,'/',num2str(N),'ingred/large.mat']);
        end
        js01 = zeros(Nsubsamps(1),Nsubsamps(2));
        js02 = zeros(Nsubsamps(1),Nsubsamps(3));
        js12 = zeros(Nsubsamps(2),Nsubsamps(3));
        for i = 1:Nsubsamps(1)
          for j = 1:Nsubsamps(2)
            js01(i,j) =...
              JSDiv(freq_large(1).subsamp_probs(:,i),maxent1_large(1).subsamp_probs(:,j));
          end
        end
        clear i j;
        for i = 1:Nsubsamps(1)
          for j = 1:Nsubsamps(3)
            js02(i,j) =...
              JSDiv(freq_large(1).subsamp_probs(:,i),maxent2_large(1).subsamp_probs(:,j));
          end
        end
        clear i j;
        for i = 1:Nsubsamps(2)
          for j = 1:Nsubsamps(3)
            js12(i,j) =...
              JSDiv(maxent1_large(1).subsamp_probs(:,i),maxent2_large(1).subsamp_probs(:,j));
          end
        end
        clear i j;
        js01_vec = reshape(js01,Nsubsamps(1)*Nsubsamps(2),1);
        js02_vec = reshape(js02,Nsubsamps(1)*Nsubsamps(3),1);
        js12_vec = reshape(js12,Nsubsamps(2)*Nsubsamps(3),1);
        JS01s(n) = mean(js01_vec);
        JS02s(n) = mean(js02_vec);
        JS12s(n) = mean(js12_vec);
        JS01_ebs(n) = std(js01_vec);
        JS02_ebs(n) = std(js02_vec);
        JS12_ebs(n) = std(js12_vec); 
      end
      clear n;
      % save the data
      Ningreds_JS = Ningreds;
      if strcmp(dataset,'cultures')
        Nrec = Nrec_used;
        save(...
          ['data/cultures/',culture_str,'/JS_',num2str(max(Ningreds)),'in.mat'],...
          'JS01s','JS02s','JS12s','JS01_ebs','JS02_ebs','JS12_ebs',...
          'Ningreds_JS','Nrec')
        allJS01s(d+k-1,:) = JS01s';
        allJS02s(d+k-1,:) = JS02s';
        allJS12s(d+k-1,:) = JS12s';
        allJS01_ebs(d+k-1,:) = JS01_ebs';
        allJS02_ebs(d+k-1,:) = JS02_ebs';
        allJS12_ebs(d+k-1,:) = JS12_ebs';
        allNrec(d+k-1) = Nrec;
        clear JS01s JS02s JS12s JS01_ebs JS02_ebs JS12_ebs;
        clear Ningred_total Ningreds_JS Nrec;
      else
        Ningred_total = Ningred;
        save(...
          ['data/',dataset,'/JS_',num2str(max(Ningreds)),'in.mat'],...
          'JS01s','JS02s','JS12s','JS01_ebs','JS02_ebs','JS12_ebs',...
          'Ningreds_JS','Nrec','Ningred_total')
        allJS01s(d+k-1,:) = JS01s';
        allJS02s(d+k-1,:) = JS02s';
        allJS12s(d+k-1,:) = JS12s';
        allJS01_ebs(d+k-1,:) = JS01_ebs';
        allJS02_ebs(d+k-1,:) = JS02_ebs';
        allJS12_ebs(d+k-1,:) = JS12_ebs';
        allNingred_total(d) = Ningred_total;
        allNrec(d+k-1) = Nrec;
        clear JS01s JS02s JS12s JS01_ebs JS02_ebs JS12_ebs;
        clear Ningred_total Ningreds_JS Nrec;
      end
    end
  end
  clear i;
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

% plot JS(P0|P1) across ingredient subsets
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(4));
shadedErrorBar(Ningreds,allJS01s(1,:)',allJS01_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
shadedErrorBar(Ningreds,allJS01s(2,:)',allJS01_ebs(2,:)',{'Color',c(2,:),'LineWidth',2});
shadedErrorBar(Ningreds,allJS01s(3,:)',allJS01_ebs(3,:)',{'Color',c(3,:),'LineWidth',2});
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(...
  [allJS01s(1,:)+allJS01_ebs(1,:),...
  allJS01s(2,:)+allJS01_ebs(2,:),...
  allJS01s(3,:)+allJS01_ebs(3,:)]);
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('D_{JS}(P^{(0)}|P^{(1)})')
text(...
  xmin+1,.95*ymax,...
  ['standard (',num2str(allNingred_total(1)),' ingredients)'],...
  'FontSize',24,'Color',c(1,:))
text(...
  xmin+1,.87*ymax,...
  ['spices (',num2str(allNingred_total(3)),' ingredients)'],...
  'FontSize',24,'Color',c(3,:))
text(...
  xmin+1,.79*ymax,...
  ['meats (',num2str(allNingred_total(2)),' ingredients)'],...
  'FontSize',24,'Color',c(2,:))
if saveflag
  file_name = ['figures/combined/JS01_insubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

% plot JS(P0|P2) across ingredient subsets
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(4));
shadedErrorBar(Ningreds,allJS02s(1,:)',allJS02_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
shadedErrorBar(Ningreds,allJS02s(2,:)',allJS02_ebs(2,:)',{'Color',c(2,:),'LineWidth',2});
shadedErrorBar(Ningreds,allJS02s(3,:)',allJS02_ebs(3,:)',{'Color',c(3,:),'LineWidth',2});
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(...
  [allJS02s(1,:)+allJS02_ebs(1,:),...
  allJS02s(2,:)+allJS02_ebs(2,:),...
  allJS02s(3,:)+allJS02_ebs(3,:)]);
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('D_{JS}(P^{(0)}|P^{(2)})')
text(...
  xmin+1,.95*ymax,...
  ['standard (',num2str(allNingred_total(1)),' ingredients)'],...
  'FontSize',24,'Color',c(1,:))
text(...
  xmin+1,.87*ymax,...
  ['spices (',num2str(allNingred_total(3)),' ingredients)'],...
  'FontSize',24,'Color',c(3,:))
text(...
  xmin+1,.79*ymax,...
  ['meats (',num2str(allNingred_total(2)),' ingredients)'],...
  'FontSize',24,'Color',c(2,:))
if saveflag
  file_name = ['figures/combined/JS02_insubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

% plot JS(P1|P2) across ingredient subsets
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(4));
shadedErrorBar(Ningreds,allJS12s(1,:)',allJS12_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
shadedErrorBar(Ningreds,allJS12s(2,:)',allJS12_ebs(2,:)',{'Color',c(2,:),'LineWidth',2});
shadedErrorBar(Ningreds,allJS12s(3,:)',allJS12_ebs(3,:)',{'Color',c(3,:),'LineWidth',2});
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(...
  [allJS12s(1,:)+allJS12_ebs(1,:),...
  allJS12s(2,:)+allJS12_ebs(2,:),...
  allJS12s(3,:)+allJS12_ebs(3,:)]);
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('D_{JS}(P^{(1)}|P^{(2)})')
text(...
  xmin+1,.95*ymax,...
  ['standard (',num2str(allNingred_total(1)),' ingredients)'],...
  'FontSize',24,'Color',c(1,:))
text(...
  xmin+1,.87*ymax,...
  ['spices (',num2str(allNingred_total(3)),' ingredients)'],...
  'FontSize',24,'Color',c(3,:))
text(...
  xmin+1,.79*ymax,...
  ['meats (',num2str(allNingred_total(2)),' ingredients)'],...
  'FontSize',24,'Color',c(2,:))
if saveflag
  file_name = ['figures/combined/JS12_insubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

% plot JS(P0|P1) across cultures
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(Ncultures+3));
shadedErrorBar(Ningreds,allJS01s(1,:)',allJS01_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
for i = 1:Ncultures
  shadedErrorBar(Ningreds,allJS01s(3+i,:)',allJS01_ebs(3+i,:)',{'Color',c(1+i,:),'LineWidth',2});
end
clear i;
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(allJS01s(1,:)+allJS01_ebs(1,:));
for i = 1:Ncultures
  ymax = max(ymax,max(allJS01s(3+i,:)+allJS01_ebs(3+i,:)));
end
clear i;
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('D_{JS}(P^{(0)}|P^{(1)})')
text(...
  xmin+1,.95*ymax,...
  ['standard (',num2str(allNrec(1)),' recipes)'],...
  'FontSize',20,'Color',c(1,:))
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
    xmin+1,(.95-.07*i)*ymax,...
    [title_str,' (',num2str(allNrec(1+i)),' recipes)'],...
    'FontSize',20,'Color',c(1+i,:))
end
clear i c title_str;
if saveflag
  file_name = ['figures/combined/JS01_recsubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

% plot JS(P0|P2) across cultures
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(Ncultures+3));
shadedErrorBar(Ningreds,allJS02s(1,:)',allJS02_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
for i = 1:Ncultures
  shadedErrorBar(Ningreds,allJS02s(3+i,:)',allJS02_ebs(3+i,:)',{'Color',c(1+i,:),'LineWidth',2});
end
clear i;
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(allJS02s(1,:)+allJS02_ebs(1,:));
for i = 1:Ncultures
  ymax = max(ymax,max(allJS02s(3+i,:)+allJS02_ebs(3+i,:)));
end
clear i;
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('D_{JS}(P^{(0)}|P^{(2)})')
text(...
  xmin+1,.95*ymax,...
  ['standard (',num2str(allNrec(1)),' recipes)'],...
  'FontSize',20,'Color',c(1,:))
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
    xmin+1,(.95-.07*i)*ymax,...
    [title_str,' (',num2str(allNrec(1+i)),' recipes)'],...
    'FontSize',20,'Color',c(1+i,:))
end
clear i c title_str;
if saveflag
  file_name = ['figures/combined/JS02_recsubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

% plot JS(P1|P2) across cultures
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(Ncultures+3));
shadedErrorBar(Ningreds,allJS12s(1,:)',allJS12_ebs(1,:)',{'Color',c(1,:),'LineWidth',2});
for i = 1:Ncultures
  shadedErrorBar(Ningreds,allJS12s(3+i,:)',allJS12_ebs(3+i,:)',{'Color',c(1+i,:),'LineWidth',2});
end
clear i;
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max(allJS12s(1,:)+allJS12_ebs(1,:));
for i = 1:Ncultures
  ymax = max(ymax,max(allJS12s(3+i,:)+allJS12_ebs(3+i,:)));
end
clear i;
ymax = 1.1*ymax;
% ymax = .35;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('D_{JS}(P^{(1)}|P^{(2)})')
text(...
  xmin+1,.95*ymax,...
  ['standard (',num2str(allNrec(1)),' recipes)'],...
  'FontSize',20,'Color',c(1,:))
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
    xmin+1,(.95-.07*i)*ymax,...
    [title_str,' (',num2str(allNrec(1+i)),' recipes)'],...
    'FontSize',20,'Color',c(1+i,:))
end
clear i culture title_str;
if saveflag
  file_name = ['figures/combined/JS12_recsubsets_',num2str(max(Ningreds)),'in'];
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

end

