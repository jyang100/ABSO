function archList2 = refinedArch(archList, A, rho)
[m,n] = size(archList);
idxDel =[];
if isempty(archList)
    archList2 = A;
    return;
end
for idx = 1:m
    p1 = archList(idx,1:n-2);
    p2 = A(1, 1:n-2);
    d = norm(p1-p2);
    if d <= rho
        if A(1,n) < archList(idx,n-1)
            A = archList(idx,:);
        end
        idxDel =[idxDel; idx];
    end
end
archList (idxDel,:)=[];
archList2 = [archList; A];
