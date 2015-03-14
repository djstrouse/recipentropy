function [] = CompareDistRE(dataset,Ningred,saveflag,titleflag,...
  showJS,showKL,meanflag,covarflag,varflag,culture)
% Scatterplot comparison of maxent model to data for top Ningred
% ingredients in dataset
%
% Written by: DJ Strouse
% Last updated: Nov 17, 2013 by DJ Strouse
%
% INPUTS
% dataset [=] string in {'standard','meats','spices','cultures'}
% Ningred [=] positive integer = number of ingredients to use
% 
% (optional)
% saveflag [=] binary = indicates whether or not to save figure
% titleflag [=] binary = indicates whether or not to title the figures
% showJS/KL [=] binary = indicates whether to display JS/KL div on plots
% mean/covar/varflag [=] binary = indicates whether to plot means/covars/vars
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
if ~exist('saveflag','var')
  saveflag = false; % don't save figs by default
end
if ~exist('titleflag','var')
  titleflag = true; % title figs by default
end
if ~exist('showJS','var') % show JS by default
  showJS = 1;
end
if ~exist('showKL','var') % show KL by default
  showKL = 1;
end
if ~exist('meanflag','var') % no means, covars, or vars by default
  meanflag = 0;
end
if ~exist('covarflag','var')
  covarflag = 0;
end
if ~exist('varflag','var')
  varflag = 0;
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
  title_str = [title_str,' (',num2str(Nrec),' recipes)'];
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
end

% load the data
if strcmp(dataset,'cultures')
  load(['data/cultures/',culture_str,'/',num2str(Ningred),'ingred/large.mat']);
else
  load(['data/',dataset,'/',num2str(Ningred),'ingred/large.mat']);
end

% set file name and prefix
if strcmp(dataset,'cultures')
  file_prefix = ['figures/cultures/',culture_str,'/'];
else
  file_prefix = ['figures/',dataset,'/']; 
end

% compare raw to maxent1
file_name = ['P0vP1_',num2str(Ningred),'in'];
compare01 = CompareDists(...
  freq_large,maxent1_large,...
  'P^{(0)}','P^{(1)}',...
  file_name,...
  'false',showJS,showKL,meanflag,covarflag,varflag);
if titleflag
  title(title_str)
end
if saveflag
  if titleflag
    file_name = [file_name,'_titled'];
  end
  saveas(compare01.hprobs,[file_prefix,file_name,'.fig']);
  export_fig([file_prefix,file_name],'-pdf','-transparent',compare01.hprobs);
end

% compare raw to maxent2
file_name = ['P0vP2_',num2str(Ningred),'in'];
compare02 = CompareDists(...
  freq_large,maxent2_large,...
  'P^{(0)}','P^{(2)}',...
  file_name,...
  'false',showJS,showKL,meanflag,covarflag,varflag);
if titleflag
  title(title_str)
end
if saveflag
  if titleflag
    file_name = [file_name,'_titled'];
  end
  saveas(compare02.hprobs,[file_prefix,file_name,'.fig']);
  export_fig([file_prefix,file_name],'-pdf','-transparent',compare02.hprobs);
end

% compare maxent1 to maxent2
file_name = ['P1vP2_',num2str(Ningred),'in'];
compare12 = CompareDists(...
  maxent1_large,maxent2_large,...
  'P^{(1)}','P^{(2)}',...
  file_name,...
  'false',showJS,showKL,meanflag,covarflag,varflag);
if titleflag
  title(title_str)
end
if saveflag
  if titleflag
    file_name = [file_name,'_titled'];
  end
  saveas(compare12.hprobs,[file_prefix,file_name,'.fig']);
  export_fig([file_prefix,file_name],'-pdf','-transparent',compare12.hprobs);
end

end

