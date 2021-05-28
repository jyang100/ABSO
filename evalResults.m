clear; clc; close all;
accuracy = 0.00001;
nr = 30;
for func_num =1:20
    fprintf('Evaluate Function %d \n', func_num);
    count = zeros(1,nr);
    ns = 0;
    for runtime = 1:nr
        % DO NOT FORGET
        initial_flag = 0; % should set the flag to 0 for each run, each function
        D = get_dimension(func_num);
        name = sprintf('results/fun_%d_run_%d.mat', func_num, runtime);
        load(name);
        % Randomize population within optimization bounds
        % (here dummy initialization within [0,1] ONLY for DEMO)
        if isempty(solutionList)
            continue;
        end
        pop = solutionList(:, 1:D);
        % How many global optima have been found?
        
        [count(runtime), goptima_found] = count_goptima(pop, func_num, accuracy);
        if count(runtime) == get_no_goptima(func_num)
            ns = ns +1;
        end
%         % Print some stuff :-)
%         if count ~=0
%             goptima_found;
%             for i=1:size(goptima_found,1)
%                 val = niching_func(goptima_found(i,:), func_num);
%                 fprintf('F_p: %f, F_g:%f, diff: %f\n', val, get_fgoptima(func_num), abs(val - get_fgoptima(func_num)))
%                 fprintf('F_p - F_g <= %f : %d\n', accuracy, abs(val - get_fgoptima(func_num))<accuracy)
%             end
%         end
    end
    pr = sum(count)/(nr*get_no_goptima(func_num));
    sr = ns/nr;
    fprintf('Peak Rate of function %d is %f \n', func_num,pr);
    fprintf('Success Rate of function %d is %f \n', func_num,sr);
    fprintf('------------------------------\n');
end