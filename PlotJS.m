function [] = PlotJS(dataset,Ningreds,saveflag,titleflag,culture)
% Plots JS divergences between P0/P1, P0/P2, and P1/P2 across Ningreds for
% dataset
%
% Written by: DJ Strouse
% Last updated: Nov 17, 2013 by DJ Strouse
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

% check if JS calculation has already been made
docalc = true;
if strcmp(dataset,'cultures')
  % first check if dir exists
  if exist(['data/cultures/',culture_str,'/JS_',num2str(max(Ningreds)),'in.mat'],'file')
    % if it does, load it
    load(['data/cultures/',culture_str,'/JS_',num2str(max(Ningreds)),'in.mat']);
    % make sure same Ningreds were used previously
    if isequal(Ningreds,Ningreds_JS)
      % if so, no need to recalc JSs
      docalc = false;
    end
  end
else
  if exist(['data/',dataset,'/JS_',num2str(max(Ningreds)),'in.mat'],'file')
    load(['data/',dataset,'/JS_',num2str(max(Ningreds)),'in.mat']);
    if isequal(Ningreds,Ningreds_JS)
      docalc = false;
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
  else
    Ningred_total = Ningred;
    save(...
      ['data/',dataset,'/JS_',num2str(max(Ningreds)),'in.mat'],...
      'JS01s','JS02s','JS12s','JS01_ebs','JS02_ebs','JS12_ebs',...
      'Ningreds_JS','Nrec','Ningred_total')
  end
end

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
    title_str = [title_str,' (',num2str(Nrec),' recipes)'];
  elseif strcmp(dataset,'meats')
    title_str = 'Meats';
    title_str = [title_str,' (',num2str(Ningred_total),' total ingredients)'];
  elseif strcmp(dataset,'spices')
    title_str = 'Spices';
    title_str = [title_str,' (',num2str(Ningred_total),' total ingredients)'];
  elseif strcmp(dataset,'standard')
    title_str = 'Standard';
    % standard compared to both ingredient and recipe subsets so label with
    % both numbers
    title_str1 = [title_str,' (',num2str(Ningred_total),' total ingredients)'];
    title_str2 = [title_str,' (',num2str(Nrec),' recipes)'];
    clear title_str;
  end
end

% plot JSs across Ningreds
xfig = 10;
yfig = 10;
wfig = 30;
hfig = 20;
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(4));
shadedErrorBar(Ningreds,JS01s,JS01_ebs,{'Color',c(1,:),'LineWidth',2});
shadedErrorBar(Ningreds,JS02s,JS02_ebs,{'Color',c(2,:),'LineWidth',2});
shadedErrorBar(Ningreds,JS12s,JS12_ebs,{'Color',c(3,:),'LineWidth',2});
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
% ymax = max([JS01s+JS01_ebs;JS02s+JS02_ebs;JS12s+JS12_ebs]);
% ymax = 1.1*ymax;
ymax = .35;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('D_{JS}')
text(xmin+1,.95*ymax,'P^{(0)} vs P^{(1)}','FontSize',24,'Color',c(1,:))
text(xmin+1,.87*ymax,'P^{(1)} vs P^{(2)}','FontSize',24,'Color',c(3,:))
text(xmin+1,.79*ymax,'P^{(0)} vs P^{(2)}','FontSize',24,'Color',c(2,:))
if titleflag
  if strcmp(dataset,'standard')
    title(title_str1)
  else
    title(title_str)
  end
end
if saveflag
  if strcmp(dataset,'cultures')
    file_name =...
      ['figures/cultures/',culture_str,'/JS_',num2str(max(Ningreds)),'in'];
  else
    file_name =...
      ['figures/',dataset,'/JS_',num2str(max(Ningreds)),'in'];
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
      ['figures/standard/JS_',num2str(max(Ningreds)),'in_titled_wNrec'];
    saveas(h,[file_name,'.fig']);
    export_fig(file_name,'-pdf','-transparent');
  end
end

end

