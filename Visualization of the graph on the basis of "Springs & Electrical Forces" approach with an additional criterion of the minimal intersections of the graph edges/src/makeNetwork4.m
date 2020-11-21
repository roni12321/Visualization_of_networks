function y = makeNetwork4 (N,K)% writes a .net+.clu files with a random network of N nodes in K clusters
    nfid = fopen ('HW2Net4.net','w');
    cfid = fopen ('HW2NetOriginal4.net','w');
    
    fprintf (nfid,'*Vertices %d\n',N); 
    fprintf (cfid,'*Vertices %d\n',N);  

    global clustersNet
    global matrixOfEdges
    
    clustVect = zeros (1,N);
    clustIndexes = zeros (1,K);
    clustersNet = zeros (2,N);
    currIndex = 1;
    matrixOfEdges = zeros (N,N,3);
    
    for i=1:N
       clustersNet(1,i) = (randi(8000)+1000)/10000;
       clustersNet(2,i) = (randi(8000)+1000)/10000;
       fprintf (cfid,'%d "ver%d" %1.4f %1.4f \n',i, i , clustersNet(1,i), clustersNet(2,i));
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
                if i~=j && randi (1000) < 9 % probability 0.9
                    lineEquation(i, j);
                    matrixOfEdges (i,j,1) = 1;
                    matrixOfEdges (j,i,1) = 1;
                end
            else
                if i~=j && randi (1000) < 10000/N % probability 10/N
                    lineEquation(i, j);
                    matrixOfEdges (i,j) = 1;
                    matrixOfEdges (j,i) = 1;
                end
            end
        end
    end
    
    % General Energy Function 
    gama1 = 0.1;
    gama2 = 0.1;
    gama3 = 0.1;
    gama4 = 0.1;
    temptemperaturehreshold = 0.01;
    clustersNetOld = zeros (2,1);
    temperature = 10;
    
    while temperature > temptemperaturehreshold
        for i=1:N
            clustersNetOld(1,1) = clustersNet(1,i);
            clustersNetOld(2,1) = clustersNet(2,i);
            
            p1 = p1Function();
            p2 = p2Function();
            p3 = p3Function();
            p4 = LineSumIntersections();
            
            uOld = gama1*p1 + gama2*p2 + gama3*p3 + gama4*p4;
            
            clustersNet(1,i) = (randi(8000)+1000)/10000;
            clustersNet(2,i) = (randi(8000)+1000)/10000;
            
            p1 = p1Function();
            p2 = p2Function();
            p3 = p3Function();
            p4 = LineSumIntersections();
            
            Unew = gama1*p1 + gama2*p2 + gama3*p3 + gama4*p4;
            
            if uOld < Unew
                p = 1-exp((uOld-Unew)/temperature); % 1-e^(U(p_old)-U(p_new)/temperatureemperature)
                if randi(1000) <= (p*1000)
                    clustersNet(1,i) = clustersNetOld(1,1);
                    clustersNet(2,i) = clustersNetOld(2,1);
                end
            end
        end
        temperature = temperature - 0.1;
    end
    
    
    for i=1:N
        fprintf (nfid,'%d "ver%d" %1.4f %1.4f 0.5000\n',i, i , clustersNet(1,i), clustersNet(2,i));
    end
    
    fprintf (nfid,'*Edges\n');
    fprintf (cfid,'*Edges\n');
    for i=1:N
        for j=i:N
            if matrixOfEdges (j,i) == 1
                fprintf (nfid,'%d %d\n',j,i);
                fprintf (cfid,'%d %d\n',j,i);
            end
        end
    end
    
    
    fclose (nfid);
    fclose (cfid);
    
    y = {clustVect',matrixOfEdges};
end


function lineEquation(i, j)

global clustersNet
global matrixOfEdges

matrixOfEdges (i,j,2) = (clustersNet(2,j)-clustersNet(2,i))/(clustersNet(1,j)-clustersNet(1,i));
matrixOfEdges (j,i,2) = matrixOfEdges (i,j,2);
matrixOfEdges (i,j,3) = clustersNet(2,i)-(matrixOfEdges (i,j,2)*clustersNet(1,i));
matrixOfEdges (j,i,3) = matrixOfEdges (i,j,3);

end

function sum = LineSumIntersections()

global matrixOfEdges
global N
sum = 0;

for i=1:N
    for j=i+1:N
        if matrixOfEdges(i,j,1) == 1
            for k=i:N
                for m=j+1:N
                   if matrixOfEdges(k,m,1) == 1 && LineIntersection(i, j, k, m)
                       sum = sum + 1;
                   end
                end
            end
        end
    end
end


end

function bool = LineIntersection(i1, j1, i2, j2)

global clustersNet
global matrixOfEdges

deltaX = matrixOfEdges (i1,j1,2) - matrixOfEdges (i2,j2,2);

if deltaX == 0
    bool = false;
else
    db = matrixOfEdges (i2,j2,3) - matrixOfEdges (i1,j1,3);
    
    bool = InternalBoundriesCheck(clustersNet(1,i1), clustersNet(1,j1), db/deltaX);
    if bool
        bool = InternalBoundriesCheck(clustersNet(1,i2), clustersNet(1,j2), db/deltaX);     
    end
    
end

end

function bool = InternalBoundriesCheck(x1, x2, ans)

bool = true;

if x1 >= x2
    if ans > x1 || ans < x2
        bool = false;
    end
else
    if ans > x2 || ans < x1
        bool = false;
    end
end

end

function sigmaV = p1Function()

global clustersNet
global N
sigmaV = 0;

for i=1:N
    for j=1:N
        if i~=j
            deltaX = (clustersNet(1,i)-clustersNet(1,j));
            deltaY = (clustersNet(2,i)-clustersNet(2,j));
            sigmaV = sigmaV + (1/(deltaX*deltaX + deltaY*deltaY));
        end
    end
end

end

function sigmaV = p2Function()

global clustersNet
global N
sigmaV = 0;

for i=1:N
    if (clustersNet(1,i)-0.1) == 0
        deltaXRight = 1;
    else
        deltaXRight = (clustersNet(1,i)-0.1);
    end
    
    if (clustersNet(1,i)-0.9) == 0
        deltaXLeft = 1;
    else
        deltaXLeft = (clustersNet(1,i)-0.9);
    end
    
    if (clustersNet(2,i)-0.1) == 0
        deltaYtop = 1;
    else
        deltaYtop = (clustersNet(1,i)-0.1);
    end
    
    if (clustersNet(2,i)-0.9) == 0
        deltaYbottom = 1;
    else
        deltaYbottom = (clustersNet(2,i)-0.9);
    end
    
    sigmaV = sigmaV + ((1/(deltaXRight*deltaXRight)) + (1/(deltaXLeft*deltaXLeft)) + (1/(deltaYtop*deltaYtop)) + (1/(deltaYbottom*deltaYbottom)));
end

end

function sigmaE = p3Function()

global clustersNet
global matrixOfEdges
global N
sigmaE = 0;

for i=1:N
    for j=1:N
        if matrixOfEdges (i,j) == 1
            deltaX = (clustersNet(1,i)-clustersNet(1,j));
            deltaY = (clustersNet(2,i)-clustersNet(2,j));
            sigmaE = sigmaE + (1/(deltaX*deltaX + deltaY*deltaY));
        end
    end
end

end