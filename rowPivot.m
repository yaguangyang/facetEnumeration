function [x,obj,iter,elapsedTime]=rowPivot(A,b,c,Ai,bi,bl,bu,rdudnt)
% A(ellxn), b(ellx1), c(1xn), Ai(mixn), bi(mix1), bl(nex1), bu(nfx1)
% If rdudnt=1, then check the redundant rows, but this cost is high
% the default is rdudnt=0
tic
if isempty(c)
    disp('Objective function must be defined');
    elapsedTime = toc;
    return;
else
    [~,n]=size(c);
    ne=n; nf=n;
end
if isempty(A)
    m=0; 
else
    [m,n]=size(A); %ne=n; nf=n;
end
if length(b)~=m
    display('rows of A and b must be the same');
    elapsedTime = toc;
    return;
end
if isempty(Ai)
    mi=0;
else
    [mi,ni]=size(Ai);
    if n~=ni
        disp('Number of column of A must be equal to the number of column of Ai'); %n=ni; ne=n; nf=n;
        elapsedTime = toc;
        return;
    end
end
if length(bi)~=mi
    display('rows of Ai and bi must be the same');
    elapsedTime = toc;
    return;
end
% Handle missing arguments
if nargin < 8
    rdudnt=0;
    if nargin < 7
        bu = [];
        if nargin < 6
            bl = [];
            if nargin < 5
                bi = [];
                if nargin < 4
                    Ai = [];
                end
            end
        end
    end
end
if isempty(bi)
    minbl=0; maxbl=0;
    minbu=0; maxbu=0;
else
    minbl=min(bl); maxbl=max(bl);
    minbu=min(bu); maxbu=max(bu);
end
if minbl==maxbl && minbl==0 && minbu==maxbu && minbu==Inf
    % standard problem
    tic;
    [x,obj,exitflag,output]=linprog(c',[],[],A,b,bl,bu);
    iter=output.iterations;
    elapsedTime=toc;
    return
end
if minbl==-Inf || minbu==Inf 
    M=10^12;
elseif minbl>-10^11 && minbu<10^11
    M=10^12;
    if minbl>-10^10 && minbu<10^10
        M=10^11;
        if minbl>-10^9 && minbu<10^9
            M=10^10;
            if minbl>-10^8 && minbu<10^8
                M=10^9;
                if minbl>-10^7 && minbu<10^7 
                    M=10^8;
                    if minbl>-10^6 && minbu<10^6 
                        M=10^7;
                        if minbl>-10^5 && minbu<10^5
                            M=10^6;
                            if minbl>-10^4 && minbu<10^4
                                M=10^5;
                                if minbl>-10^3 && minbu<10^3
                                    M=10^4;
                                    if minbl>-10^2 && minbu<10^2 
                                        M=10^3; 
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

if isempty(bl)  % the default lower bound is zero
    bl=zeros(n,1);
else
    for ii=1:n
        if bl(ii)==-Inf 
            bl(ii)=-M;
        end
    end
end
if isempty(bu)  % the default upper bound is M
    bu=M*ones(n,1);
else
    for ii=1:n    
        if bu(ii)==Inf
            bu(ii)=M;
        end
    end
end
epsilon=10^(-6);
yc=c; ell=bl; E=eye(n); F=eye(n); 
E0=[];       % index set of eq consraints in Base
if m>0
    E1=1:m;  % index set of eq consraints in Base
elseif m==0
    E1=[];
end
I0=1:n; I1=n+1:2*n+mi; %initial I0=EIdx, I1=[ne+FIdx,ne+nf+AIIdx]
LI0=length(I0);
for i=1:n
    if c(i)<0
        yc(i)=-c(i);
        E(i,i)=-1; F(i,i)=-1;
        ell(i)=-bu(i); bu(i)=-bl(i);
    end
end
F=-E; bf=-bu; % G=[E; F];
% remove redundant non-base constraints
fdel=find(bf==-M);
F(fdel,:)=[]; bf(fdel)=[];
fdel=find(bf==M);
F(fdel,:)=[]; bf(fdel)=[];
nf=length(bf); I1=n+1:n+nf+mi;
AIdx=1:m; EIdx=1:ne; FIdx=1:nf; AIIdx=1:mi;
% remove redundant non-base constraints
ba=ell;
x=E*ell;  % because inv(E)=E
iter=0;
obj=c*x; 
if isempty(b)
    sigmE=[];
else
    sigmE=A*x-b; 
end
sigmL=E*x-ell; sigmU=F*x-bf; %sigmU=F*x+bf;
if isempty(bi)
    sigmI=[];
else
    sigmI=Ai*x-bi; 
end
if m==0
    maxsigE=0;
elseif m>0
    [maxsigE,idE]=max(abs(sigmE));
end
[maxsigL,idL]=min(sigmL);
[maxsigU,idU]=min(sigmU);
if mi==0
    maxsigI=0;
elseif mi>0
    [maxsigI,idI]=min(sigmI);
end
Ab=E;

while maxsigE/max(1,norm(x))>epsilon || maxsigL<-epsilon ...
        || maxsigU<-epsilon || maxsigI/max(1,norm(x))<-epsilon
    % remove the redundant equality constraint if there is any
    reCal=0;
    if ~isempty(E1)
        for ii=length(E1):-1:1
            yd=LUsolution(Ab',A(E1(ii),:)');% yd=Ab'\A(AIdx(ii),:)';
            if abs(A(E1(ii),:)*x-b(E1(ii)))<=10^(-10) ...
                    && max(abs(yd(1:LI0)))<=10^(-10) 
                E1(ii)=[]; sigmE(ii)=[];
                AIdx(ii)=[];
                if ii<idE
                    idE=idE-1;
                elseif ii==idE
                    reCal=1;
                end
                m=m-1;
            end
        end
    end
    % recalculate maxsigE
    if reCal==1 
        if isempty(E1)
            maxsigE=0; idE=0;
        else
            [maxsigE,idE]=max(abs(sigmE));
        end
    end
    if rdudnt==1
        % remove the redundant inequality constraint in E if there is any
        reCal=0;
        for ii=ne:-1:1
            if ismember(EIdx(ii),I1) %&& ii~=idL
                yd=LUsolution(Ab',E(EIdx(ii),:)');% yd1=Ab'\E(EIdx(ii),:)';
                if E(EIdx(ii),:)*x>ell(EIdx(ii)) ... %-10^10 ...
                        && min(yd(1:LI0))>=0 % -10^10 
                    I1=setdiff(I1,EIdx(ii)); sigmL(ii)=[];
                    EIdx(ii)=[];
                    if ii<idL
                        idL=idL-1;
                    elseif ii==idL
                        reCal=1;
                    end
                    ne=ne-1;
                end
            end
        end
        if reCal==1 && maxsigE<10^(-10)
            if ne==0
                maxsigL=0; idL=0;
            elseif ne>0
                [maxsigL,idL]=min(sigmL);
            end
        end
        % remove the redundant inequality constraint in F if there is any
        reCal=0;
        for ii=nf:-1:1
            if ismember(FIdx(ii)+n,I1) %&& ii~=idU
                yd=LUsolution(Ab',F(FIdx(ii),:)');% yd=Ab'\F(FIdx(ii),:)';
                if F(FIdx(ii),:)*x>bf(FIdx(ii)) ... %-10^10 ...
                        && min(yd(1:LI0))>=0 %-10^10
                    I1=setdiff(I1,n+FIdx(ii)); sigmU(ii)=[];
                    FIdx(ii)=[];
                    if ii<idU
                        idU=idU-1;
                    elseif ii==idU
                        reCal=1;
                    end
                    nf=nf-1;
                end
            end
        end
        if reCal==1 && maxsigE<10^(-10)
            if nf==0
                maxsigU=0; idU=0;
            elseif nf>0
                [maxsigU,idU]=min(sigmU);
            end
        end
        % remove the redundant inequality constraint in Ai if there is any
        reCal=0;
        if mi>0
            for ii=mi:-1:1
                if ismember(AIIdx(ii)+2*n,I1) && ii~=idI
                    yd=LUsolution(Ab',Ai(AIIdx(ii),:)');%yd=Ab'\Ai(AIIdx(ii),:)';
                    if Ai(AIIdx(ii),:)*x>bi(AIIdx(ii)) ... %-10^10 ...
                            && min(yd(1:LI0))>=0 %-10^10 
                        I1=setdiff(I1,n+n+AIIdx(ii));  sigmI(ii)=[];
                        AIIdx(ii)=[];
                        if ii<idI
                            idI=idI-1;
                        elseif ii==idI
                            reCal=1;
                        end
                        mi=mi-1;
                    end
                end
            end
        end
        if reCal==1 && maxsigE<10^(-10)
            if mi==0
                maxsigI=0; iDI=0;
            elseif mi>0
                [maxsigI,idI]=min(sigmI);
            end
        end
    end
    % select entering row, virtual A matrix is in the order of equality,
    % lower and upper bounds, and inequality constraints
    if ~isempty(E1) && maxsigE>=10^(-10)
        eV=1; ap=A(E1(idE),:); bp=b(E1(idE));
    else
        [minsig,idC]=min([maxsigL maxsigU maxsigI]);
        if minsig<0 && idC==1
            eV=0; ap=E(EIdx(idL),:); bp=ell(EIdx(idL));
        elseif minsig<0 && idC==2
            eV=-1; ap=F(FIdx(idU),:); bp=bf(FIdx(idU));
        elseif minsig<0 && idC==3
            eV=-2; ap=Ai(AIIdx(idI),:); bp=bi(AIIdx(idI));
        end
    end
    yp=LUsolution(Ab',ap');  
    zidx=find(abs(yp)<10^(-5));  % 10^(-10));
    yp(zidx)=0;
    % if there is no feasible solution, exit while loop
    % idE will be calculated at the end of the program
    if ~isempty(E1) && m>0 && eV==1
        if ((ap*x)>bp+10^(-8) && min(yp(1:LI0))>=-10^(-8))     
            warning('There is no feasible solution')
            infeasible=1;
            elapsedTime = toc;
            return;
        elseif ((ap*x)<bp-10^(-8) && max(yp(1:LI0))<=10^(-8))      
            warning('There is no feasible solution')
            infeasible=1;
            elapsedTime = toc;
            return;
        end
    end
    % if there is no feasible solution, exit while loop
    % idL, idU, and idI constraints will be considered later
    if ~isempty(I0) % && ne>0  && eV==0
        if ((ap*x)<bp-10^(-8) && max(yp(1:LI0))<=0)
            warning('There is no feasible solution')
            infeasible=1;
            elapsedTime = toc;
            return;
        end
    end
    % select leaving row and record base and nonbase rows
    if eV==1  && m>0      % entering variable p in E1 equality constraints
        if ap*x<bp && max(yp(1:LI0))>0 % Case 1
        tmpInitial=M;
            for ii=1:LI0
                if yp(ii)>0
                    tmp=yc(ii)/yp(ii);
                    if tmp<tmpInitial
                        q=ii;
                        tmpInitial=tmp;
                    end
                end
            end
        elseif ap*x>bp && min(yp(1:LI0))<0 % Case 2
            tmpInitial=-M;
            for ii=1:LI0
                if yp(ii)<0
                    tmp=yc(ii)/yp(ii);
                    if tmp>tmpInitial
                        q=ii;
                        tmpInitial=tmp;
                    end
                end
            end
        end
        E0=[E0 E1(idE)]; %LE0=length(E0);   % E0 are equality constraint in base, equality E1(idE) enters base
        E1(idE)=[]; % E1=setdiff(E1,E1(idE));  
%         I1=[I1 I0(q)]; remove artificial leaving constraint
        if I0(q)<=n && ell(I0(q))==-M
            EIdx=setdiff(EIdx,I0(q)); ne=ne-1;
        else
            I1=[I1 I0(q)];
        end
        I0(q)=[];  LI0=LI0-1;   % I0 are inequality constraints in base. Inequality constraint q leaves base 
        ba(q)=[];
        ba(n)=b(AIdx(idE)); %%%ba(n)=b(E1(idE));
        Ab(q,:)=[]; Ab=[Ab; ap]; % Ab=[E(I0,:); A(E0,:)];
    elseif eV==0 % entering variable p from matrix E % I=(lower bound; upper bound; inequality)
        tmpInitial=M;
        for ii=1:LI0
            if yp(ii)>0
                tmp=yc(ii)/yp(ii);
                if tmp<tmpInitial
                    q=ii;
                    tmpInitial=tmp;
                end
            end
        end
        I1=setdiff(I1,EIdx(idL));     % leave an inequality constraint
        I1=[I0(q) I1];
        I0(q)=[]; 
        I0=[EIdx(idL) I0];            % enter an inequality constraint
        ba(q)=[];
        ba=[ell(EIdx(idL)); ba];
        Ab(q,:)=[];
        Ab=[ap; Ab]; % Ab=[E(EIdx(idL),:); Ab];
    elseif eV==-1 % entering variable p from matrix F % I=(lower bound; upper bound; inequality)
        tmpInitial=M;
        for ii=1:LI0
            if yp(ii)>0
                tmp=yc(ii)/yp(ii);
                if tmp<tmpInitial
                    q=ii;
                    tmpInitial=tmp;
                end
            end
        end
        I1=setdiff(I1,FIdx(idU)+n);     % leave an inequality constraint
        I1=[I0(q) I1];
        I0(q)=[]; 
        I0=[FIdx(idU)+ne I0];            % enter an inequality constraint
        % LI0=length(I0);         the length does not change
        ba(q)=[];
        ba=[bf(FIdx(idU)); ba];
        Ab(q,:)=[];
        Ab=[ap; Ab]; % Ab=[F(FIdx(idU),:); Ab];
    elseif eV==-2  && mi>0 % entering variable p from matrix Ai % I=(lower bound; upper bound; inequality)
        tmpInitial=10^6;
        for ii=1:LI0
            if yp(ii)>0
                tmp=yc(ii)/yp(ii);
                if tmp<tmpInitial
                    q=ii;
                    tmpInitial=tmp;
                end
            end
        end
        I1=setdiff(I1,AIIdx(idI)+2*n);     % leave an inequality constraint
        I1=[I0(q) I1];
        I0(q)=[]; 
        I0=[AIIdx(idI)+n+nf I0];            % enter an inequality constraint
        ba(q)=[];
        ba=[bi(AIIdx(idI)); ba];
        Ab(q,:)=[];
        Ab=[ap; Ab]; % Ab=[Ai(AIIdx(idI),:); Ab];
    end
    % Compute c^(k+1)
    coef=yc(q)/yp(q);
    if abs(coef)<10^(-10)
        coef=0;
    end
    for ii=1:n
        if ii~=q
            yc(ii)=yc(ii)-yp(ii)*coef;
            if abs(yc(ii))<10^(-10)
                yc(ii)=0;
            end
        end
    end
    yc(q)=[];
    if eV==1 % equality entering row is place in the bolttom of Ab
        yc(n)=coef;
    else     % inequality entering row is place in the top of Ab
        yc=[coef yc];
    end
    % Remove redundant leaving row
    % not implemented because netlib problems don't have  redundant ineq
    % compute x variable
    x=LUsolution(Ab,ba);
    % Compute the constraint violation and obj ...
    if isempty(E1)
        sigmE=[];
    else
        sigmE=A(E1,:)*x-b(E1); 
    end
    if isempty(EIdx)
        sigmL=[];
    else
        sigmL=E(EIdx,:)*x-ell(EIdx); 
    end
    if isempty(FIdx)
        sigmU=[];
    else
        sigmU=F(FIdx,:)*x-bf(FIdx);
    end
    if isempty(bi)
        sigmI=[];
    else
        sigmI=Ai(AIIdx,:)*x-bi(AIIdx); 
    end
    if m==0 || length(E0)==m
        maxsigE=0; idE=0;
    elseif m>0
        [maxsigE,idE]=max(abs(sigmE));
    end
    if ne==0
        maxsigL=0; idL=0;
    elseif ne>0
        [maxsigL,idL]=min(sigmL);
    end
    if nf==0
        maxsigU=0; idU=0;
    elseif nf>0
        [maxsigU,idU]=min(sigmU);
    end
    if mi==0
        maxsigI=0; iDI=0;
    elseif mi>0
        [maxsigI,idI]=min(sigmI);
    end
    iter=iter+1;    
    if mod(iter,100)==0 
        disp(['iteration: ', num2str(iter)]);
    end
    obj=yc*ba; % obj=c*x; formula (25)
end

for i=1:n
    if x(i)>-epsilon && x(i)<0
        x(i)=0;
    end
end
if max(abs(ba))==M
    display('the solution is unbounded');
    for i=1:n
        if ba(i)==M
            ba(i)=inf;
        elseif ba(i)==-M
            ba(i)=-inf;
        end
    end
    obj=yc*ba;
end
% toc
elapsedTime = toc;

function [x]=LUsolution(Ab,ba)
epsilon=10^(-5);
[L,U]=lu(Ab); 
zrow=find(abs(diag(U))<epsilon);
ti=length(zrow);
sbar=L\ba;
if ti>0
    for kk=1:ti
        U(zrow(kk),zrow(kk))=epsilon; 
    end
end
x=U\sbar;
    