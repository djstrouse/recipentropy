close all; clear all;
cd('/Users/djstrouse/Documents/Dropbox/recipentropy');
for n = [2 3]
  for o = 0:2
    for s = 1:2
      cd('/Users/djstrouse/Documents/Dropbox/recipentropy');
      RecipEntropySingle(n,o,s);
    end
    clear s;
  end
  clear o;
end
clear n;
cd('/Users/djstrouse/Documents/Dropbox/recipentropy');