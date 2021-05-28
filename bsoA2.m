function archList= bsoA2(func_num,population_size,number_dimension,range_left,range_right,max_iteration, rho)

%probability_dimension_disruption = 0.2; % probability to determine whether a dimension is disrupted or not
probability_disrupt = 0.1; % probability for disrupting elitists. one elitis every 5 generations, and only one dimension;
probability_elitist = 0.8; % probability for select elitist, not normals, to generate new individual;
probability_one = 0.8; % probability for select one individual, not two, to generate new individual;
percentage_elitist = 0.1;
logsig_slope = 25; % slope of the s-shape function


fitness = ones(population_size,max_iteration);
attn = zeros(population_size,max_iteration);
pos = cell(population_size,number_dimension,max_iteration);


reDisNum =5;
%archList = zeros(1,number_dimension+1);

if number_dimension == 1
    x = range_left:0.001:range_right;
    y = -1*ones(1,size(x,2));
    for idx = 1:size(x,2)
        y(idx) =  niching_func(x(idx), func_num);
    end
    plot(x,y, 'linewidth',2);
%     xlabel("x_1")
%     ylabel("Fitness value");
    hold on;
end

if number_dimension == 2
    x1 = range_left(1):0.01:range_right(1);
    x2 = range_left(2):0.01:range_right(2);
    y = -1*ones(size(x1,2), size(x2,2));
    for idx1 = 1:size(x1,2)
        for idx2 = 1:size(x2,2)
            y(idx1,idx2) = niching_func([x1(idx1) x2(idx2)],func_num);
        end
    end
    popPlot = surf(x1,x2,y');
    axis tight;
    set(popPlot,'LineStyle','none');
%     xlabel("x_1")
%     ylabel("x_2")
%     zlabel("Fitness value")
    colorbar;
    hold on;
end

% restart = 1;
archList = double.empty;
archList2 = double.empty;

current_iteration_number = 1; % initialize current iterantion  number
% restartIter  = 1;
archListCount = 1;
% start the main loop of the BSO algorithm
max_iteration_comp = ceil(max_iteration - (2*number_dimension+1)*(max_iteration * probability_disrupt)/population_size); % compensate the fitness evaulation due to disruption
%max_iteration_comp = max_iteration - (max_iteration/population_size); % compensate the fitness evaulation due to disruption
% most_atten = ones(max_iteration_comp,1);  % store best fitness value for each iteration

population = range_left + (range_right - range_left) .* rand(population_size,number_dimension); % initialize the populationlation of individuals
population_temporary  = range_left + (range_right - range_left) .* rand(population_size,number_dimension); % initialize the temporary populationlation of individuals
number_elitist = round(population_size * percentage_elitist); % number of elitists
number_normal = population_size - number_elitist;             % number of normals
%  index_original = ones(population_size,1); % initialize corresponding original indexs in the population of sorted fitness values
% best_fitness = ones(max_iteration,1);  % store best fitness value for each iteration
fitness_population = ones(population_size,1);  % store fitness value for each individual in each population
attn_population = zeros(population_size,1);  % store fitness value for each individual in each population
% fitness_population_sorted = ones(population_size,1);  % store  sorted fitness values
% attn_population_sorted = ones(population_size,1);  % store  sorted fitness values
individual_temporary = zeros(1,number_dimension);  % store a temporary individual
% calculate fitness for each individual in the initialized populationlation
for idx = 1:population_size
    fitness_population(idx,1) = niching_func(population(idx,:),func_num);
%     attn_population(idx,1) = getAttn(func_num, population(idx,:), fitness_population(idx,1),rho);
end

% popTotal = population;
% fitTotal = fitness_population;

for idx = 1:population_size
%     [attn_population(idx,1),~] = getNearAttn(population(idx,:), fitness_population(idx,1), popTotal, fitTotal);
    [attn_population(idx,1),~] = getNearAttn(population(idx,:), fitness_population(idx,1), population, fitness_population, attn_population, archList);
end

while current_iteration_number <= max_iteration_comp % use the compesented maximum iteration
%     step_size = rho*ones(1,number_dimension);
    step_size =  logsig(((0.5*max_iteration_comp - current_iteration_number)/logsig_slope)) * rand(1,number_dimension);
    %     if restart == 1
    %         step_size = 0.1*ones(1,number_dimension); % effecting the step size of generating new individuals by adding random values
    
    
    %         restartIter = 0;
    %         step_size = logsig(((0.5*max_iteration - restartIter)/logsig_slope)) * rand(1,number_dimension);
    %         restart = 0;
    %     end
    %     eps = step_size;
    % disrupt individual
    if current_iteration_number > 1 % don't do disruption for the first generation
        %disrupt every genration but for one dimension of one individual
        r_1 = rand();  % generate a randon value
        if r_1 < probability_disrupt % decide whether to select one individual to be disrupted
            idx = ceil(population_size * rand()); % index of the selected individual
            individual_temporary(1,:) = population(idx,1,:); % temporary individual = selected individual
            selected_dimension = ceil(number_dimension * rand());
            individual_temporary(1,selected_dimension) = range_left(selected_dimension) + (range_right(selected_dimension) - range_left(selected_dimension)) * rand(); %one dimention of selected individual to be replaceed by a random number
            fv = niching_func(individual_temporary,func_num); % evaluate the disrupted individual
%             [av, flag] = getAttn(func_num, individual_temporary, fv,rho);

            [av, flag] = getNearAttn(individual_temporary, fv, population, fitness_population, attn_population, archList);
%             [av, flag] = getNearAttn(population(idx,:), fitness_population(idx,1), popTotal, fitTotal);
            
            population(idx,:) = individual_temporary(1,:); % replace the selected individual with the disrupted one
            
            fitness_population(idx,1) = fv; % assign the fitness value to the disrupted individual
            attn_population(idx,1) = av; % assign the attention value to the disrupted individual
            
%             popTotal = [popTotal; individual_temporary(1,:)];
%             fitTotal = [fitTotal;fv];
            
            if flag == 1 % if most salient
                A = [individual_temporary(1,:) av fv];
                archList2 = refinedArch(archList, A, rho);
            end
        end
    end

    if archListCount > reDisNum
        population(index_original(number_elitist+1:population_size),:) = range_left + (range_right - range_left) .* rand(number_normal,number_dimension);
        for idx = number_elitist+1:population_size
            fitness_population(idx,1) = niching_func(population(idx,:),func_num);
            %             [attn_population(idx,1), flag] = getAttn(func_num, individual_temporary, fv,rho);
            %             [attn_population(idx,1),flag] = getNearAttn(individual_temporary, fv, population, fitness_population);
%             popTotal = [popTotal; population(idx,:)];
%             fitTotal = [fitTotal;fitness_population(idx,1)];
        end
        
        for idx = number_elitist+1:population_size
            [attn_population(idx,1),~] = getNearAttn(population(idx,:), fitness_population(idx,1), population, fitness_population, attn_population, archList);
%              [attn_population(idx,1),~] = getNearAttn(population(idx,:), fitness_population(idx,1), popTotal, fitTotal);
        end
%         restartIter  = 1;
        archListCount = 1;
    end
    
      % sort individuals in a population based on their attention values
    [attn_population_sorted, index_original] = sort(attn_population,'descend');
%     record the best fitness in each iteration
%     most_atten(current_iteration_number, 1) = attn_population_sorted(1,1);
    
    % generate population_size new individuals by adding Gaussian random values
    for idx = 1:population_size
        % form the seed individual
        r_1 = rand();  % generate a randon value
        if r_1 < probability_elitist % select elitists to generate a new individual
%             if size(archList,1) > 2 && rand() < 0.5
%                 r = rand();
%                 if r < probability_one % use one elitist to generate a new individual
%                     archIdx = randsample(size(archList,1), 1);
%                     individual_temporary(1,:) = archList(archIdx,1:number_dimension);
%                 else % use two elitists to generate a new individual
%                     tem = rand();
%                     archIdx = randsample(size(archList,1), 2);
%                     individual_temporary(1,:) = tem * archList(archIdx(1),1:number_dimension) + (1-tem) * archList(archIdx(2),1:number_dimension);
%                 end
                %             if rand() < 0.5 && ~isempty(archList)
                %                 archIdx = randsample(size(archList,1), 1);
                %                 individual_temporary(1,:) = archList(archIdx,1:number_dimension);
                %                 step_size = logsig(((0.5*max_iteration - current_iteration_number)/logsig_slope)) * rand(1,number_dimension);
%             else
                r = rand(); % generate a random number
                inx_selected_one = ceil(number_elitist * rand());
                inx_selected_two = ceil(number_elitist * rand());
                if r < probability_one % use one elitist to generate a new individual
                    individual_temporary(1,:) = population(index_original(inx_selected_one,1),:);
                else % use two elitists to generate a new individual
                    tem = rand();
                    individual_temporary(1,:) = tem * population(index_original(inx_selected_one,1),:) + (1-tem) * population(index_original(inx_selected_two,1),:);
                end
%             end
        else  % select normals to generate a new individual
            r = rand(); % generate a random number
            inx_selected_one = number_elitist + ceil(number_normal * rand());
            inx_selected_two = number_elitist + ceil(number_normal * rand());
            if r < probability_one % use one normal to generate a new individual
                individual_temporary(1,:) = population(index_original(inx_selected_one,1),:);
            else % use two elitists to generate a new individual
                tem = rand();
                individual_temporary(1,:) = tem * population(index_original(inx_selected_one,1),:) + (1-tem) * population(index_original(inx_selected_two,1),:);
            end
        end
        
        % add Gaussian random value to seed individual to generate a new individual
        
        a = individual_temporary(1,:) + step_size .* normrnd(0,1,1,number_dimension);
        while min(a(1,:) - range_left) < 0 || max(a(1,:)- range_right) > 0
            a = range_left + (range_right - range_left) .* rand(1,number_dimension);
%             a = individual_temporary(1,:) + step_size .* normrnd(0,1,1,number_dimension);
        end
        individual_temporary(1,:) = a;
        
        % selection between new one and the old one with the same index
        fv = niching_func(individual_temporary,func_num);
%         [av, flag] = getAttn(func_num, individual_temporary, fv,rho);
                 [av, flag] = getNearAttn(individual_temporary, fv, population, fitness_population, attn_population,archList);
%                  [av, flag] = getNearAttn(population(idx,:), fitness_population(idx,1), popTotal, fitTotal);
                 
%                  popTotal = [popTotal; individual_temporary(1,:)];
%                  fitTotal = [fitTotal;fv];
       
% %        if  av > attn_population(idx,1)  % better than the previous one, replace
%        if  fv > fitness_population(idx,1)   % better than the previous one, replace
%            attn_population(idx,1) = av;
%            fitness_population(idx,1) = fv;
%            population_temporary(idx,:) = individual_temporary(1,:);
%        else % keep the previous one
%            population_temporary(idx,:) = population(idx,:);
%        end
       
        nearIdx = knnsearch(population, individual_temporary, 'k', 1);
        if  av > attn_population(nearIdx,1)  % better than the nearest one, replace
%          if  fv > fitness_population(nearIdx,1)   % better than the nearest one, replace
            attn_population(nearIdx,1) = av;
            fitness_population(nearIdx,1) = fv;
            population_temporary(nearIdx,:) = individual_temporary(1,:);
            
        else % keep the previous one
            population_temporary(nearIdx,:) = population(nearIdx,:);
        end
        if flag == 1
            A = [individual_temporary(1,:) av fv];
            archList2 = refinedArch(archList, A, rho);
%             archList = archList2;
        end
    end
    
    % copy temporary population to population to start next iteration
    for idx = 1:population_size
        population(idx,:) = population_temporary(idx,:);
    end
    
    if size(archList2,1) == size(archList,1)
         archListCount = archListCount +1;
%         if archList2 == archList
%             archListCount = archListCount +1;
%         else
%             archListCount = 1;
%         end
    else
        archListCount = 1;
    end
%     archList = archList2;
    archList = refineArch(archList2,rho);
    
    
    % update main loop counter
    fitness(:,current_iteration_number) = fitness_population;
    attn(:,current_iteration_number) = attn_population;
    pos(:,:,current_iteration_number) = num2cell(population);
    
    current_iteration_number = current_iteration_number +1;
%     restartIter = restartIter + 1;
    
end

if isempty(archList)
    archList = [population,fitness_population];
end
% archList = [population, fitness_population];
% filename = 'test.mat';
% save(filename);
if number_dimension == 1 || number_dimension == 2
    %         if current_iteration_number >1
    %             delete(popPlot);
    %         end
    if ~isempty(archList)
        archPos = archList(:,1:(size(archList,2)-2));
        archFit = archList(:,size(archList,2));
        if number_dimension == 1
            scatter(archPos, archFit,'filled', 'red');
            name = sprintf('results/f%dr.eps', func_num);
            saveas(gcf,name,'psc2');
%             xlabel("x_1")
%             ylabel("Fitness value");
        else
            scatter3(archPos(:,1), archPos(:,2), archFit,'filled', 'red');
            colorbar;
%             xlabel("x_1")
%             ylabel("x_2");
            view(0,90);
             name = sprintf('results/f%dr.eps', func_num);
             saveas(gcf,name,'psc2');
        end
    end
    
    
    maxFit =  zeros(1,max_iteration_comp);
    maxAttn =   zeros*ones(1,max_iteration_comp);

    for i = 1:max_iteration_comp
        if i==1
            maxFit(i:max_iteration_comp) = max(fitness(:,1));
            maxAttn(i:max_iteration_comp) = max(attn(:,1));
            continue;
        end
        ma = max(attn(:,i));
        ma1 = max(maxAttn);
        mf = max(fitness(:,i));
        mf1 =  max(maxFit);
        if mf < mf1
            maxFit(i:max_iteration_comp) = mf1;
        else
            maxFit(i:max_iteration_comp) = mf;
        end
        
        if ma < ma1
            maxAttn(i:max_iteration_comp) = ma1;
        else
            maxAttn(i:max_iteration_comp) = ma;
        end
    end
    figure(2)
    plot(maxFit,'linewidth',2);
    name = sprintf('results/f%drf.eps', func_num);
    saveas(gcf,name,'psc2');
    figure(3)
    plot(maxAttn,'linewidth',2);
    name = sprintf('results/f%dra.eps', func_num);
    saveas(gcf,name,'psc2');
    close all;
end
end



