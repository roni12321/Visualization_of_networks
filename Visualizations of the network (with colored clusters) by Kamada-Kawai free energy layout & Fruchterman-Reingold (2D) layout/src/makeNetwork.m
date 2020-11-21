function y = makeNetwork (N,K)% writes a .net+.clu files with a random network of N nodes in K clusters
    nfid = fopen ('HW2Net.net','w');
    cfid = fopen ('HW2Clust.clu','w');
    
    fprintf (nfid,'*Vertices %d\n',N); 
    fprintf (cfid,'*Vertices %d\n',N);  

    clustVect = zeros (1,N);
    clustIndexes = zeros (1,K);
    currIndex = 1;
    edgeMatrix = zeros (N,N);
    
    for i=1:(K-1)
        currIndex = currIndex + N/K + randi (ceil((N*0.5)/K)) - (N*0.25)/K; % base cluster size is N/K +- 0.25 N/K
        clustIndexes(1,i) = currIndex;
    end
    clustIndexes(1,K) = N;
    
    currIndex = 1;
    for i=1:N
        clustVect (1,i) = currIndex;
        fprintf (cfid,'%d\n',currIndex);
        if (i > clustIndexes(1,currIndex))
            currIndex = currIndex + 1;
        end
    end
    
    fprintf (nfid,'*Edges\n'); 
    for i=1:N
        for j=i:N
            if clustVect(1,i) == clustVect(1,j) 
                if i~=j && randi (1000) < 900 % probability 0.9
                     edgeMatrix (i,j) = 1;
                     edgeMatrix (j,i) = 1;
                end
            else
                if i~=j && randi (1000) < 10000/N % probability 10/N
                    edgeMatrix (i,j) = 1;
                    edgeMatrix (j,i) = 1;
                end
            end
            if edgeMatrix (j,i) == 1
                fprintf (nfid,'%d %d\n',j,i); 
            end
        end
    end
    
    fclose (nfid);
    fclose (cfid);
    
    y = {clustVect',edgeMatrix};
end