function [] = PlotKL(dataset,Ningreds,saveflag,titleflag,culture)
% Plots KL divergences between P0/P1, P0/P2, and P1/P2 across Ningreds for
% dataset
%
% NEEDS UPDATED TO JS STANDARDS BEFORE USE
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
  elseif strcmp(dataset,'meats')
    title_str = 'Meats';
  elseif strcmp(dataset,'spices')
    title_str = 'Spices';
  elseif strcmp(dataset,'standard')
    title_str = 'Standard';
  end
  title_str = [title_str,' (',num2str(Nrec),' recipes)'];
end

% load and process the data
KL01s = zeros(length(Ningreds),1);
KL01_ebs = zeros(length(Ningreds),1);
KL02s = zeros(length(Ningreds),1);
KL02_ebs = zeros(length(Ningreds),1);
KL12s = zeros(length(Ningreds),1);
KL12_ebs = zeros(length(Ningreds),1);
for n = 1:length(Ningreds)
  N = Ningreds(n);
  if strcmp(dataset,'cultures')
    load(['data/cultures/',culture_str,'/',num2str(N),'ingred/large.mat']);
  else
    load(['data/',dataset,'/',num2str(N),'ingred/large.mat']);
  end
  kl01 = zeros(Nsubsamps(1),Nsubsamps(2));
  kl02 = zeros(Nsubsamps(1),Nsubsamps(3));
  kl12 = zeros(Nsubsamps(2),Nsubsamps(3));
  for i = 1:Nsubsamps(1)
    for j = 1:Nsubsamps(2)
      kl01(i,j) =...
        KLDiv(freq_large(1).subsamp_probs(:,i),maxent1_large(1).subsamp_probs(:,j));
    end
  end
  clear i j;
  for i = 1:Nsubsamps(1)
    for j = 1:Nsubsamps(3)
      kl02(i,j) =...
        KLDiv(freq_large(1).subsamp_probs(:,i),maxent2_large(1).subsamp_probs(:,j));
    end
  end
  clear i j;
  for i = 1:Nsubsamps(2)
    for j = 1:Nsubsamps(3)
      kl12(i,j) =...
        KLDiv(maxent1_large(1).subsamp_probs(:,i),maxent2_large(1).subsamp_probs(:,j));
    end
  end
  clear i j;
  kl01_vec = reshape(kl01,Nsubsamps(1)*Nsubsamps(2),1);
  kl02_vec = reshape(kl02,Nsubsamps(1)*Nsubsamps(3),1);
  kl12_vec = reshape(kl12,Nsubsamps(2)*Nsubsamps(3),1);
  KL01s(n) = mean(kl01_vec);
  KL02s(n) = mean(kl02_vec);
  KL12s(n) = mean(kl12_vec);
  KL01_ebs(n) = std(kl01_vec);
  KL02_ebs(n) = std(kl02_vec);
  KL12_ebs(n) = std(kl12_vec); 
end
clear n;

% plot KLs across Ningreds
xfig = 10;
yfig = 10;
wfig = 30;
hfig = 20;
h = figure;
hold off
hold on
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
c = colormap(pmkmp(4));
shadedErrorBar(Ningreds,KL01s,KL01_ebs,{'Color',c(1,:),'LineWidth',2});
shadedErrorBar(Ningreds,KL02s,KL02_ebs,{'Color',c(2,:),'LineWidth',2});
shadedErrorBar(Ningreds,KL12s,KL12_ebs,{'Color',c(3,:),'LineWidth',2});
prettyplot(24)
xmin = min(Ningreds)-1;
xmax = max(Ningreds)+1;
xlim([xmin xmax])
ymax = max([KL01s+KL01_ebs;KL02s+KL02_ebs;KL12s+KL12_ebs]);
ymax = 1.1*ymax;
ylim([0 ymax])
set(gca,'XTick',Ningreds)
xlabel('number of ingredients')
ylabel('D_{KL}')
if titleflag
  title(title_str)
end
text(xmin+1,.95*ymax,'P^{(0)} vs P^{(1)}','FontSize',24,'Color',c(1,:))
text(xmin+1,.87*ymax,'P^{(0)} vs P^{(2)}','FontSize',24,'Color',c(2,:))
text(xmin+1,.79*ymax,'P^{(1)} vs P^{(2)}','FontSize',24,'Color',c(3,:))
if saveflag
  if strcmp(dataset,'cultures')
    file_name =...
      ['figures/cultures/',culture_str,'/KL_',...
      num2str(max(Ningreds)),'in'];
  else
    file_name =...
      ['figures/',dataset,'/KL_',...
      num2str(max(Ningreds)),'in'];
  end
  if titleflag
    file_name = [file_name,'_titled'];
  end
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

end

