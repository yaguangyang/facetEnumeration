%              Yaguang Yang 
%              May 22, 2025
%  Copyright (c) 2025 by Y, Yang, all rights reserved.
%
% This function performs facet enumeration. It implements Algorithm of
% [Yang2025]. User will provide vertex matrix (represent vertex in row),
% and the function will return facet expressions in inequality constraints.
% Several examples are included. Users can uncommment any example and
% run the program to get the result.
% 
% Users may change the code, but the copywight will still belong to Y. Yang.
% If any report/publication is based on this code, the authors should 
% acknowledge the use and cite [Yang2025].
% 
% References:
% [Yang2025] Yaguang Yang, A Facet Enumeration Algorithm for Convex
% Polytopes, arXiv:1909.11843, 2025. 
% https://doi.org/10.48550/arXiv.1909.11843

clear all
epsilon = 10^(-8); Bgst=10^(10);

% the first problem is a triangle, test good
A = [ 0 0; 3 0; 0 3]; 
% % the second problem is a cube, test good
% A=[ 0 0 0;
%    0 0 1;
%    0 1 0;
%    1 0 0;
%    0 1 1;
%    1 0 1;
%    1 1 0;
%    1 1 1];
% % the third problem is a octahedon, test good
%  A=[ 0  0  1;
%     -1  0  0;
%      0 -1  0;
%      1  0  0;
%      0  1  0;
%      0  0  -1];
% % the fourth problem is a cross-polytope, test good
%  A=[ 0  0  0  1;
%     -1  0  0  0;
%      0 -1  0  0;
%      0  0 -1  0;
%      1  0  0  0;
%      0  1  0  0;
%      0  0  1  0;
%      0  0  0  -1];
% % the fifth problem is given by a reviewer
% A=[1 0 1;
%     1 0 -1;
%     1.25 -1 1;
%     1.25 -1 -1;
%     0.25 -1 1;
%     0.25 -1 -1;
%     -1 0 1;
%     -1 0 -1;
%     -1.25 1 1;
%     -1.25 1 -1;
%     -0.25 1  1;
%     -0.25 1 -1];
% % the sixth problem is a simplex in n-dimensional space, max(n)<=9
% n=5;
% A0=zeros(1,n); A1=(n+1)*eye(n);
% A=[A0; A1];
% % the seventh problem is a regular icosahedron (20 facets\12 vertices\30 edges)
% % https://www.wolframalpha.com/input/?i=regular+icosahedron+vertex+coordinates
% a=sqrt(50-10*sqrt(5)); b=sqrt(2/(5-sqrt(5))); c=sqrt(10-2*sqrt(5));
% d=1+sqrt(5); dm=1-sqrt(5); ep=5+sqrt(5); em=5-sqrt(5);
% A=[0 0 -5/a;
%    0 0 5/a;
%    -b 0 -1/c;
%    b 0 1/c;
%    d/(2*c) -1/2 -1/c;
%    d/(2*c) 1/2 -1/c;
%    -d/(2*c) -1/2 1/c;
%    -d/(2*c) 1/2 1/c;
%    dm/(2*c) -0.5*sqrt(ep/em) -1/c;
%    dm/(2*c) 0.5*sqrt(ep/em) -1/c;
%    -dm/(2*c) -0.5*sqrt(ep/em) 1/c;
%    -dm/(2*c) 0.5*sqrt(ep/em) 1/c]
% % the eighth problem has 3n facets, 2n+2 vertices, and 5n edges.
% n=8; p=360/n;
% deg=p*pi/180;
% A1=zeros(n,3); A2=zeros(n,3); A=[0 0 -2];
% for ii=1:n
%     A1(ii,:)=[cos((ii-1)*deg), sin((ii-1)*deg), -1];
%     A2(ii,:)=[cos((ii-1)*deg), sin((ii-1)*deg), 1];
% end
% A=[A; A1; A2]; 
% A=[A; [0 0 2]];

% starting timer
tic
% compute centroid
meanA=mean(A);
% move polytope center to origin
[n,d] = size(A);
Ac=zeros(n,d); % Am=Ac; % Ac is shifted polytope
e=ones(d,1);
warning('off','all');
for ii=1:n
    Ac(ii,:)=A(ii,:)-meanA;
end

nHplane=0;
% find the neighbor vertices for every vertex
D=zeros(n,n); % ith row saves neighbor vertices of the ith vertex
% F=[];         % ith row saves vertices on the ith facet 
for ii=1:n
    for jj=ii+1:n
        ujmi=Ac(jj,:)-Ac(ii,:);
        t=-Ac(ii,:)*ujmi'/(norm(ujmi))^2;
        z=Ac(ii,:)+t*ujmi;
        normZ=norm(z);
        if normZ >epsilon  %%% z~=0 therefore it is not at the origin
            uijd2=(Ac(ii,:)+Ac(jj,:))/2;
            Ain=Ac;
            Ain([ii,jj],:)=[];
            normAi=zeros(n-2,1);
            for kk=1:n-2
                Ain(kk,:)=Ain(kk,:)-uijd2; % u_k-(u_i+u_j)/2
                normAi(kk)=norm(Ain(kk,:));
            end
            bin=-Ain*uijd2';                % (u_k-(u_i+u_j)/2)*((u_i+u_j)/2)'
            Ain=-[Ain, -normAi];            % [u_k-(u_i+u_j)/2, - \| . \|]
            cin=[zeros(d,1); 1];
            Aeq=[ujmi 0];                  % [u_j - u_i 0]
            beq=uijd2*ujmi';               % ((u_i+u_j)/2)(u_j - u_i)'
            lb=[-ones(d,1)+uijd2'; -10^10];
            ub=[ones(d,1)+uijd2'; 10];
            if ii==1 && jj==6
                [ii jj]
            end
%             [yf,fval,exflag,oput] = linprog(cin,Ain,bin,Aeq,beq,lb,ub);
            [yf,fval,iter,elapsedTime]=rowPivot(Aeq,beq,cin',Ain,bin,lb,ub);

            if fval<-epsilon
                D(ii,jj)=1;
                D(jj,ii)=1;
            else
                D(ii,jj)=0;
                if abs(fval)<=epsilon && norm(yf(1:3)-uijd2')> epsilon %remark 4

                elseif abs(fval)<=epsilon && norm(yf(1:3)-uijd2')<= epsilon
                    % if obj=0 and y=r  % remark 2
                    Ain=[Ain; e' 0];
                    bin=[bin; 1+sum(uijd2)];
                    [yf,fval,iter,elapsedTime]=rowPivot(Aeq,beq,cin',Ain,bin,lb,ub);
                    if fval<-epsilon
                        D(ii,jj)=1;
                        D(jj,ii)=1;
                    end
                end
            end
        else
            %%%%%{ii,jj} cross the oridge therefore is not an edge
        end
    end
end

H=[];         % total number of hyperplanes found
B=[];         % corresponding branched
Si=zeros(d,d); 
Blength=0; examinedPath=0;
toc

U0=[]; Uc=1; Ut=2:n;     % initialize U_0 and U_t
while length(U0) < n
    Uc2=[];              % Uc2 save vertices which will be in Uc in next iteration
    for ii=1:length(Uc)  % for each vertex in Uc, do the following
%         display(Uc)
%         display(ii)
        levelRow=zeros(d,1);
        Done=0;
        Si=zeros(d,d); 
        Si(1,:)=Ac(Uc(ii),:);
        lastRow=Uc(ii);
        D(:,lastRow)=selectV(D(:,lastRow)); % color the column not usable (replace all 1 by 2 in column lastRow of D)
        level=1;
        levelRow(level)=Uc(ii);  % save the vertex of this level
        while Done==0;
            fOne=firstMinInRow(D(lastRow,:)); % find the next adjacent vertex (the first index fOne of D(lastRow,fOne)=1)
            if level<d & fOne>0 % if the length of the tree under Uc(ii)<d
                level=level+1;  % the length of the tree + 1
                Uc2=[Uc2 fOne];  % the newly found adjacent vertex should be in next Uc
                Si(level,:)=Ac(fOne,:);  % the base matrix adds the vertex
                levelRow(level)=fOne;    % the tree member is updated
                D(lastRow,fOne)=2;
                lastRow=fOne;
                D(:,lastRow)=selectV(D(:,lastRow)); %color vertices in D (replace element=1 in column lastRow by 2)
            elseif level==d
                sortL=sort(levelRow);
                newPath=1;
                for jj=1:Blength
                    if sortL==B(:,jj)
                        newPath=0;
                        break;
                    end
                end
                if rcond(Si)>epsilon && newPath==1
                    examinedPath=examinedPath+1;
                    h=Si\e;
                    rept=0;                    % 0 =new hyperplane is found
                    if ~isempty(H)
                        for kk=1:nHplane
                            if norm(H(:,kk)-h)<epsilon
                                rept = 1;      % 1= the hyperplane was found
                                break;
                            end
                        end
                    end
                    if rept==0    % 0 =new hyperplane is found
                        fisbility = length(find(Ac*h > 1+epsilon));  % 0=feasible, 1=infeasible
                        corct = Ac(Uc(ii),:)*h;           % corct>=1 valid, corct<1 invalid
%                         F1=zeros(1,n);
                        if fisbility==0 && corct>=1-epsilon        % corct>=1?
                            H=[H h];
                            B=[B sortL];   %B may be used to check if new h is in B
                            Blength=Blength+1;
                            nHplane=nHplane+1;
                        end
                    end
                end
                levelRow(level)=0;  % move up one level in the vertex tree
                level=level-1;
                lastRow=levelRow(level);
                fIdx=find(D(lastRow,:)); % get index in lastRow row which is positive
                fMin=min(D(lastRow,fIdx)); % min in fOne row and fIdx colums

                while fMin==2  % if there is no more unexamined vertex, 
                    % move up one level and re-color the adjacent matrix
                    if level>1
                        D(levelRow(level-1),levelRow(level))=3; % set item = 3, 
                        % item = 3, all vertex below the node has been searched 
                        levelRow(level)=0;
                        level=level-1;
                        lastRow=levelRow(level);
                        D=resetD(D); 
                        for kk=1:level
                           D(:,levelRow(kk))=selectV(D(:,levelRow(kk))); % replace 1 by 2
                           % in column vector levelRow(kk) 
                        end
                        fIdx=find(D(levelRow(level),:)); % get index in fOne row
                        fMin=min(D(lastRow,fIdx)); % min in fOne row and fIdx colums
                    else
                        D=resetD(D);
                        fMin=min(D(lastRow,fIdx));
                        Done=1;
                    end
                    if fMin==2
                        D(levelRow(level),:)=d3to1(D(levelRow(level),:));
                    elseif fMin==3 % every vertex in this row has been checked
                        D(lastRow,:)=d3to1(D(lastRow,:));
                        Done=1;
                    end
                end

            elseif level==0 & fOne==0
                Done=1;
            end
        end
        Uc2=removeDuplicate(Uc2); 
        D=resetD(D);
    end % for Uc(ii)

    U0=[U0 Uc];
    U0=removeDuplicate(U0);
    Uc2=removeDuplicate(Uc2); 
    idx = ismember(Uc2,U0);
    Uc2(idx) = [];
    idx = ismember(Uc2,Uc);
    Uc2(idx) = [];
    Uc=Uc2;
    clear Uc2;
end   % while
%--------
% Line 14 shift the polytope
%--------
b=zeros(size(H,2),1);
for ii=1:size(H,2)    %recover the center at original place
    b(ii)=1+meanA*H(:,ii);
end
[H' b]
timer=toc