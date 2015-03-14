function [] = zipf(dataset,Ningred,logtype,normflag,saveflag,titleflag,culture)
% Makes a zipf plot using Ningred most frequent ingredients from dataset
%
% Written by: DJ Strouse
% Last updated: Nov 16, 2013 by DJ Strouse
%
% INPUTS
% dataset [=] string in {'standard','meats','spices','cultures'}
% Ningred_used [=] positive integer = number of ingredients to use
% 
% (optional)
% logtype [=] string in {'linear','semilog','loglog'}
% normflag [=] binary = indicates whether yscale should show prob/freq
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
% add option for placing ingredient names as xlabels

% init
if ~exist('logtype','var')
  logtype = 'linear';
end
if ~exist('normflag','var')
  normflag = false;
end
if ~exist('saveflag','var')
  saveflag = false;
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

% load the data
if strcmp(dataset,'standard')
  load('data/recipes.mat'); clear dict recipes recipes_binary;
elseif strcmp(dataset,'meats')
  load('data/meats.mat'); clear meat_dict final_meat_binary;
  ingred_freq = meat_ingred_freq; clear meat_ingred_freq;
  ingreds = final_meat_ingred; clear final_meat_ingred;
elseif strcmp(dataset,'spices')
  load('data/spices.mat'); clear spice_dict final_spices_binary;
  ingred_freq = ingred_freq_spices; clear ingred_freq_spices;
  ingreds = final_spices; clear final_spices;
elseif strcmp(dataset,'cultures')
  load('data/recipes.mat'); clear dict ingred_freq;
  [~,recipes_binary]=...
  ExtractCulturalRecipes(culture,recipes,recipes_binary); clear recipes;
  ingred_freq = sum(recipes_binary,1); clear recipes_binary;
else
  error('dataset string not recognized!')
end
[sorted_ingred_freqs,i] = sort(ingred_freq,'descend');
sorted_ingred_freqs = sorted_ingred_freqs(1:Ningred);
to_use = i(1:Ningred); clear i;
sorted_ingred_names = ingreds(to_use); clear to_use ingreds;

% plot
xfig = 10;
yfig = 10;
wfig = 30;
hfig = 20;
buffer = .05;

if normflag
  yvar = sorted_ingred_freqs./sum(sorted_ingred_freqs);  
  normflagstr = 'normflag';
else
  yvar = sorted_ingred_freqs;
  normflagstr = 'counts';
end

h = figure;
set(h,'units','centimeters','outerposition',[xfig yfig wfig hfig])
if strcmp(logtype,'linear')
  plot(1:Ningred,yvar,'.','LineWidth',2);
elseif strcmp(logtype,'semilog')
  semilogy(1:Ningred,yvar,'.','LineWidth',2);
elseif strcmp(logtype,'loglog')
  loglog(1:Ningred,yvar,'.','LineWidth',2);
else
  error('logtype not a recognizable string!')
end
% set(gca,'XTick',1:Ningred)
prettyplot(20)
% xlim([0 Ningred+1])
% ylim([yvar(Ningred)-buffer*yvar(Ningred),yvar(1)+buffer*yvar(1)])
xlabel('ingredient frequency rank')
ylabel('ingredient frequency')
if titleflag
  title(title_str)
end
if saveflag
  if strcmp(dataset,'cultures')
    file_name =...
      ['figures/cultures/',culture_str,'/zipf_',...
      num2str(Ningred),'in',logtype,'_',normflagstr];
  else
    file_name =...
      ['figures/',dataset,'/zipf_',...
      num2str(Ningred),'in',logtype,'_',normflagstr];
  end
  if titleflag
    file_name = [file_name,'_titled'];
  end
  saveas(h,[file_name,'.fig']);
  export_fig(file_name,'-pdf','-transparent');
end

end

