function [nearAttn, flag]= getNearAttn(ind, fv, pop, popFit, popAttn, archList)
if ~isempty(archList)
    [~,n] = size(archList);
    pop = [pop;archList(:,1:n-2)];
    popAttn = [popAttn;archList(:,n-1)];
    popFit = [popFit;archList(:,n)];
end
    nearIdx = knnsearch(pop, ind, 'k', 20);
%       nearIdx = 1:100;
      nearAttn = 0;
      flag = 1;
      for i = 1:length(nearIdx)
          a = popFit(nearIdx(i));
          nearAttn = nearAttn + (fv-a)*exp(-norm(ind - pop(nearIdx(i),:))/0.4); 
          if flag ==1 && (fv < a)
              flag = 0;
          end
      end
      nearAttn = nearAttn/length(nearIdx);
      
%        for i = 1:length(nearIdx)
%            if popAttn(nearIdx(i)) == 0
%                flag = 0;
%                return;
%            end
%            if flag ==1 && nearAttn <= popAttn(nearIdx(i))
%                flag = 0;
%                return;
%            end
%        end
end