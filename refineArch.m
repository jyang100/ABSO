function archList = refineArch(archList, rho)
[m,n] = size(archList);
idxDel =[];
if isempty(archList)
    return;
end
for idx = 1:m
    p1 = archList(idx,1:n-2);
    for idx2= 1:m
        if idx2 == idx
            continue;
        end
    p2 = archList(idx2,1:n-2);
    d = norm(p1-p2);
    if d <= rho
        if archList(idx2,n) < archList(idx,n)
            idxDel =[idxDel; idx2];
        else
            idxDel =[idxDel; idx];
        end
    end
    end
end
archList (idxDel,:)=[];
