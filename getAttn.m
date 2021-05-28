function [av, flag] = getAttn(func_num, x, fv, rho)
% rho = 0.1;
n = length(x);
flag = 1;
 sr_1 = 0;
 sr_2 = 0;
for idx = 1:n
    x1 = x;
    x2 = x;
    xeps1 = x(idx) + rho;
    xeps2 = x(idx) - rho;
    x1(idx) = xeps1;
    x2(idx) = xeps2;
    fv1 = niching_func(x1,func_num);
    fv2 = niching_func(x2,func_num);
    if flag ==1 && (fv < fv1 || fv < fv2)
            flag = 0;
    end
    sr_1 =  sr_1 + fv1;
    sr_2 =  sr_2 + fv2;
end
sr = sr_1 + sr_2;
av = 2*n*fv - sr;
end