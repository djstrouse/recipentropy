%check that files exist

types = {'meats','spices','standard'};
orders = [0 1 2];
subsamps = 1:25;
n_ingreds = [4 6 8 10 12 14 16 18];
dir = '/Users/arianasp/Desktop/Dropbox/recipentropy/della_runs_output';

disp('List of files that do not exist:')
for type_idx = 1:length(types)
    type = types{type_idx};
    for n_ingred = n_ingreds
        for order = orders
            for subsamp = subsamps
                name = [dir '/' type '/' num2str(n_ingred) 'ingred/order' num2str(order) '/subsamp' num2str(subsamp) '/subsample.mat'];
                if not(exist(name))
                    disp(name)
                end
            end
        end
    end
end

types = {'EastAsian','LatinAmerican','NorthAmerican','SouthernEuropean','WesternEuropean'};
orders = [0 1 2];
subsamps = 1:25;
n_ingreds = [4 6 8 10 12 14 16 18];
dir = '/Users/arianasp/Desktop/Dropbox/recipentropy/della_runs_output/cultures';

disp('List of files that do not exist:')
for type_idx = 1:length(types)
    type = types{type_idx};
    for n_ingred = n_ingreds
        for order = orders
            for subsamp = subsamps
                name = [dir '/' type '/' num2str(n_ingred) 'ingred/order' num2str(order) '/subsamp' num2str(subsamp) '/subsample.mat'];
                if not(exist(name))
                    disp(name)
                end
            end
        end
    end
end
        