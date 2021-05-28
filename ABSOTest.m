% ******************************************************************************
% * Version: 1.0.1
% * Last modified on: 4 April, 2016 
% * Developers: Michael G. Epitropakis, Xiaodong Li.
% *      email: mge_(AT)_cs_(DOT)_stir_(DOT)_ac_(DOT)_uk 
% *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
% * ****************************************************************************

% Demonstration file on how to use the benchmark suite of the Competition 
% Bellow you can find examples on how to: 
%       i) evaluate a solution, and (niching_func)
%      ii) calculate the number of global optima in a set of solutions (count_goptima)
% 
clear; clc; close all;

%DO NOT FORGET
global initial_flag; % the global flag used in test suite 
NP=100;
accuracy = 0.1;
nr = 1;
for func_num =1:1
    for runtime = 1:nr
        name = sprintf('results/fun_%d_run_%d.mat', func_num, runtime);
%         if exist(name,'file')
%             continue;
%         end
        fprintf('running_%d of fun_%d \n', runtime, func_num);
        % Set the lower and upper bound for each function
        % DO NOT FORGET
        initial_flag = 0; % should set the flag to 0 for each run, each function
        D = get_dimension(func_num);
        MaxFes = get_maxfes(func_num);
%         max_iteration  = ceil((MaxFes/NP)/(2*D+1));
        max_iteration  = MaxFes/NP;
        % Dimension of the problem
        
        range_left = get_lb(func_num);
        range_right=get_ub(func_num);
        rho = get_rho(func_num);
        
        % Potential solution
        solutionList = bsoA2(func_num,NP,D,range_left,range_right,max_iteration, rho);
        
        % Save the solution
        
%         parsave(name,solutionList);
        
        % 	val = niching_func(x, func_num); % fitness evaluation
        % 	fprintf('f_%d : f(1...1) = %f\n',func_num,solutionList);
    end
end


% % 
% for func_num =1:1
% 	% DO NOT FORGET
% 	initial_flag = 0; % should set the flag to 0 for each run, each function 
% 	D = get_dimension(func_num);
% 	MaxFes = get_maxfes(func_num);
% 
% 	% Randomize population within optimization bounds 
% 	% (here dummy initialization within [0,1] ONLY for DEMO)
% % 	pop = rand(NP,D);
%     pop = solutionList(:, 1:D);
% 	% How many global optima have been found?
% 
% 	[count, goptima_found] = count_goptima(pop, func_num, accuracy);
% 	fprintf('f_%d, In the current population there are %d global optima!\n', func_num,count);
%     get_no_goptima(func_num)
% 
% 	% Print some stuff :-)
% 	if count ~=0
% 		goptima_found;
% 		for i=1:size(goptima_found,1)
% 			val = niching_func(goptima_found(i,:), func_num);
% 			fprintf('F_p: %f, F_g:%f, diff: %f\n', val, get_fgoptima(func_num), abs(val - get_fgoptima(func_num)))
% 			fprintf('F_p - F_g <= %f : %d\n', accuracy, abs(val - get_fgoptima(func_num))<accuracy)
% 		end
% 	end
% end
