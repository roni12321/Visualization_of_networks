function y = makeNetwork3 (N,K)% writes a .net+.clu files with a random network of N nodes in K clusters
    nfid = fopen ('HW2Net.net','w');
    cfid = fopen ('HW2NetOriginal.net','w');
    
    fprintf (nfid,'*Vertices %d\n',N); 
    fprintf (cfid,'*Vertices %d\n',N);  

    clustVect = zeros (1,N);
    clustIndexes = zeros (1,K);
    clustNet = zeros (2,N);
    currIndex = 1;
    edgeMatrix = zeros (N,N);
    
    for i=1:N
       clustNet(1,i) = (randi(8000)+1000)/10000;
       clustNet(2,i) = (randi(8000)+1000)/10000;
       fprintf (cfid,'%d "ver%d" %1.4f %1.4f 0.5000\n',i, i , clustNet(1,i), clustNet(2,i));
    end
    
    for i=1:(K-1)
        currIndex = currIndex + N/K + randi (ceil((N*0.5)/K)) - (N*0.25)/K; % base cluster size is N/K +- 0.25 N/K
        clustIndexes(1,i) = currIndex;
    end
    clustIndexes(1,K) = N;
    
    currIndex = 1;
    for i=1:N
        clustVect (1,i) = currIndex;
        %fprintf (cfid,'%d\n',currIndex);
        if (i > clustIndexes(1,currIndex))
            currIndex = currIndex + 1;
        end
    end
    
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
        end
    end
    
    k = 1;
    Delta = 0.4;
    
    Const = 1.6 * (10^-19);
    Colon = 8.98755 * (10^9);
    
    Flag = 1;
    
    while Flag == 1
        Flag = 0;
        for i=1:N

            IfMovFlag = 1;
            while IfMovFlag == 1
                IfMovFlag=0;
                SumFx = 0;
                SumFy = 0;
            
                for j=1:N
%Electric-force
                    if i~=j
                        DeltaX = (clustNet(1,i)-clustNet(1,j));
                        DeltaY = (clustNet(2,i)-clustNet(2,j));
                        Func = (Colon*(Const*Const))/(DeltaX*DeltaX + DeltaY*DeltaY);
                        Angle = atand(abs(DeltaY/DeltaX));

                        SumFx = SumFx + cosd(Angle)*Func*DeltaX;
                        SumFy = SumFy + sind(Angle)*Func*DeltaY;
                    end
                        
%Spring-force
                    if edgeMatrix (i,j) == 1

                        DeltaX = (clustNet(1,i)-clustNet(1,j));
                        DeltaY = (clustNet(2,i)-clustNet(2,j));
                        Func = ((-1) * k)*(Delta-sqrt(DeltaX*DeltaX + DeltaY*DeltaY));
                        Angle = atand(abs(DeltaY/DeltaX));
                        
                        if DeltaX >= 0
                            DeltaX = -1;
                        else
                            DeltaX = 1;
                        end
                        
                        if DeltaY >= 0
                            DeltaY = -1;
                        else
                            DeltaY = 1;
                        end
                        
                        SumFx = SumFx + cosd(Angle)*Func*DeltaX;
                        SumFy = SumFy + sind(Angle)*Func*DeltaY;
                        
                        
                    end
                end

                if SumFx > 0.01
                    if SumFx > 0.02
                        Flag = 1;
                        IfMovFlag=1;
                    end
                    clustNet(1,i) = clustNet(1,i) + 0.001; 
                end
                if SumFx < -0.01
                    if SumFx < -0.02
                        Flag = 1;
                        IfMovFlag=1;
                    end
                    clustNet(1,i) = clustNet(1,i) - 0.001;
                end
                if SumFy > 0.01
                    if SumFy > 0.02
                        Flag = 1;
                        IfMovFlag=1;
                    end
                   clustNet(2,i) = clustNet(2,i) + 0.001;
                end
                if SumFy < -0.01
                    if SumFy < -0.02
                        Flag = 1;
                        IfMovFlag=1;
                    end
                    clustNet(2,i) = clustNet(2,i) - 0.001;
                end
                
				if abs(SumFx) + abs(SumFy) <= 0.02
					Flag = 0
					IfMovFlag=0
				end
            end
        end
    end
    
    for i=1:N
        fprintf (nfid,'%d "ver%d" %1.4f %1.4f 0.5000\n',i, i , clustNet(1,i), clustNet(2,i));
    end
    
    fprintf (nfid,'*Edges\n');
    fprintf (cfid,'*Edges\n');
    for i=1:N
        for j=i:N
            if edgeMatrix (j,i) == 1
                fprintf (nfid,'%d %d\n',j,i);
                fprintf (cfid,'%d %d\n',j,i);
            end
        end
    end
    
    
    fclose (nfid);
    fclose (cfid);
    
    y = {clustVect',edgeMatrix};
end