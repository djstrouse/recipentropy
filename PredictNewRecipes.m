%% init

close all; clear all;
load('data/recipentropy_16in_large');
load('data/recipes.mat');

%% find Ningred_used most common ingredients

Ningred_used = log2(length(maxent1_large.probs));
[sorted_ingred_freq,i] = sort(ingred_freq,'descend');
to_use = i(1:Ningred_used); clear i;
recipes_binary_subsamp = recipes_binary(:,to_use)';
ingreds_subsamp = ingreds(to_use); clear to_use;
disp(sprintf('Using the following %i ingredients:',Ningred_used));
disp(ingreds_subsamp');

%% make sorted list of most probable antirecipes

antirecipe_indices = find(ingred_freq_subsamp_large.probs==0);
antirecipe_predprobs = maxent2_large(1).probs(antirecipe_indices);
Narecs = length(antirecipe_predprobs);
[sorted_antirecipe_predprobs,antirecipe_sorting_indices] = sort(antirecipe_predprobs,'descend');

%% display k most probable antirecipes

k = 5;
predicted_recs_binary = zeros(k,16);
for j = 1:k
  predicted_recs_binary(j,:) =...
    int2bin(antirecipe_indices(antirecipe_sorting_indices(j))-1,16);
  prob = sorted_antirecipe_predprobs(j);
  eb = maxent2_large(1).prob_ebs(antirecipe_indices(antirecipe_sorting_indices(j)));
  [decplace,str] = sigdig(eb);
  disp(...
    [sprintf(['The %ith most probable (p = ',str],j,prob)...
    char(177),...
    sprintf([str,') antirecipe is:'],eb)])
  disp(ingreds_subsamp(find(predicted_recs_binary(j,:)))')
  clear prob eb;
end
clear j k;

%% display k least probable antirecipes

k = 5;
predicted_recs_binary = zeros(k,16);
for j = 1:k
  index = Narecs-j+1;
  predicted_recs_binary(index,:) =...
    int2bin(antirecipe_indices(antirecipe_sorting_indices(index))-1,16);
  prob = sorted_antirecipe_predprobs(index);
  eb = maxent2_large(1).prob_ebs(antirecipe_indices(antirecipe_sorting_indices(index)));
  [decplace,str] = sigdig(eb);
  disp(...
    [sprintf(['The %ith least probable (p = ',str],j,prob)...
    char(177),...
    sprintf([str,') antirecipe is:'],eb)])
  disp(ingreds_subsamp(find(predicted_recs_binary(index,:)))')
  clear prob eb index;
end
clear j;

%% make sorted list of most underrepresented recipes

recipe_indices = find(sum(ingred_freq_subsamp_large.subsamp_probs~=0,2)==25);
Nrecs = length(recipe_indices);
Nsubsamp = ingred_freq_subsamp_large.Nsubsamp;
freqs = ingred_freq_subsamp_large.subsamp_probs(recipe_indices,:);
probs = maxent2_large(1).subsamp_probs(recipe_indices,:);
subsamp_ratios = zeros(Nrecs,Nsubsamp);
for j = 1:Nsubsamp
  subsamp_ratios(:,j) = probs(:,j)./freqs(:,j);
end
clear j;
ratios_mean = mean(subsamp_ratios,2);
ratios_std = std(subsamp_ratios,0,2);
[sorted_ratios recipe_sorting_indices] = sort(ratios_mean,'descend');

%% display k most underrepresented recipes

k = 5;
predicted_recs_binary = zeros(k,16);
for j = 1:k
  predicted_recs_binary(j,:) =...
    int2bin(recipe_indices(recipe_sorting_indices(j))-1,16);
  ratio = sorted_ratios(j);
  eb = ratios_std(recipe_sorting_indices(j));
  [decplace,str] = sigdig(eb);
  disp(...
    [sprintf(['The %ith most underrepresented (prob/freq = ',str],j,ratio)...
    char(177),...
    sprintf([str,') recipe is:'],eb)])
  disp(ingreds_subsamp(find(predicted_recs_binary(j,:)))')
  clear ratio eb;
end
clear j;

%% display k most overrepresented recipes

k = 5;
predicted_recs_binary = zeros(k,16);
for j = 1:k
  index = Nrecs-j+1;
  predicted_recs_binary(index,:) =...
    int2bin(recipe_indices(recipe_sorting_indices(index))-1,16);
  ratio = sorted_ratios(index);
  eb = ratios_std(recipe_sorting_indices(index));
  [decplace,str] = sigdig(eb);
  disp(...
    [sprintf(['The %ith most overrepresented (prob/freq = ',str],j,ratio)...
    char(177),...
    sprintf([str,') recipe is:'],eb)])
  disp(ingreds_subsamp(find(predicted_recs_binary(index,:)))')
  clear ratio eb index;
end
clear j;