clear;
clc;
N = 357;
net = makeNetwork (N,7);
a = net{2};


fileName = 'HW2Clust_ex.clu';

netT = a; % temporary network matrix, contains only rows/cols related to our currently selected cluster
lapla = netT;

aa = sum (a);
for i=1:length(netT)
	lapla(i,i)= sum (netT(:,i)); %build laplacian matrix
end
    
lapla = lapla - netT * 2; % make all weight negative
[egVect,egDig] = eig (lapla);
clustVectT = egVect(:,2); % take lambda2 vect
for j=1:length(clustVectT)
    if (clustVectT(j) < 0) % modify the partition depending on lambda2
        clustVect(j) = 1;
    else
        clustVect(j) = 2;
    end
end

cfid = fopen (fileName,'w');
fprintf (cfid,'*Vertices %d\n',N);  
    
for i=1:N
	fprintf (cfid,'%d\n',clustVect(i));
end

fclose (cfid);
