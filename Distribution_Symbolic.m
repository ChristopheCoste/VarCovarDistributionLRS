clear all; clc;close all;

%% inputs
s=2; % number of states

S=[.3 0; .2 .7] % expectation of survival

x = sym('x_%d',[1 s]); % dummy variable of length s
y = sym('y'); %dummy univariate variable

pgfF=[x(1)/2 + 1/2, x(1)/2 + x(2)/2] % pgf of reproduction 

%% intermediary outputs

pgfS  = ones(1,s) +  (x- ones(1,s))*S; %pgf of survival (eq 4)
vpa(pgfS,2)


%F=[0.5 0.5; 0 0.5] %expectation of reproduction


F=expect_pgf(pgfF,x,s,s) % expectation of reproduction
VF=Var_pgf(pgfF,x,s,s)% variance covaraince of reproduction
% these two elements are not necessary to compute the pgf of LRS

%% pgf n1
pgfn0=x(1)*x(2)

pgfn1=subs(pgfn0,x,pgfF.*pgfS)

expand(pgfn1)

%% pgf of R

D=(eye(s)-S*diag(pgfF))
D^(-1)

pgfR=(pgfF + ones(1,s)*D - ones(1,s))*D^(-1) %eq 14

%% pgf of total LRS (rT)
pgff=subs(pgfF,x,y*ones(1,s)); %pgf of univariate (i.e. unstructured) reproduction
dd=diag(pgff);
pgfr=(pgff-ones(1,s)*S*dd)*(eye(s)-S*dd)^(-1); 


%% Distributions by succesive derivations of the pgfs

maxn=15;


% distributions for unstructured LRS : rT
distres1=dist_pgf(pgfr(1),y,maxn); %distribution of rT (small type)
distres2=dist_pgf(pgfr(2),y,maxn); %distribution of rT (medium type)

% distributions for structured LRS 
Dist=zeros(maxn+1,maxn+1,s); %distribution for R
for i=1:s
Der=pgfR(i);
Derx1=Der;
for n1=0:maxn
Derx2=Derx1;
for n2=0:maxn
Dist(n1+1,n2+1,i)=(1/factorial(n1))*(1/factorial(n2))*subs(Derx2,{x(1),x(2)},{0,0});
Derx2=diff(Derx2,x(2));
end
Derx1=diff(Derx1,x(1));
end
end

%% figure

nmax3d=8

figure
t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
nexttile;
bar3(Dist(1:nmax3d+1,1:nmax3d+1,1),'b')
xlabel('type \emph{small} offspring','Interpreter','latex')
ylabel('type \emph{large} offspring','Interpreter','latex')
zlabel('probability','Interpreter','latex')
title('(a): joint distribution of LRS of \emph{small} type: ${\mathbf{R_1}}$','Interpreter','latex')
view(-135,20)
set(gca,'XTickLabel',[0:nmax3d])
set(gca,'YTickLabel',[0:nmax3d])

nexttile;
bar3(Dist(1:nmax3d+1,1:nmax3d+1,2),'r')
xlabel('type \emph{small} offspring','Interpreter','latex')
ylabel('type \emph{large} offspring','Interpreter','latex')
zlabel('probability','Interpreter','latex')
title('(b): joint distribution of LRS of \emph{large} type: ${\mathbf{R_2}}$','Interpreter','latex')
view(-135,20)
set(gca,'XTickLabel',[0:nmax3d])
set(gca,'YTickLabel',[0:nmax3d])

nexttile;
bar(distres1(1,:),distres1(2,:),'b')
xlabel('lifetime reproductive success','Interpreter','latex')
ylabel('probability','Interpreter','latex')
title('(c): distribution of total LRS of \emph{small} type: $r(1)$','Interpreter','latex')

nexttile;
bar(distres2(1,:),distres2(2,:),'r')
xlabel('lifetime reproductive success','Interpreter','latex')
ylabel('probability','Interpreter','latex')
title('(d): distribution of total LRS of \emph{large} type: $r(2)$','Interpreter','latex')

%% functions


% expectation from pgf: from first derivative valued at 1
function [z]=expect_pgf(pgf,x,q,qi)  %means the pgf is of size qi with elements in x_1,...,x_qo
%a=length(pgf);
expect_pgf=zeros(q,qi);
for i=1:q
expect_pgf(i,:)=(subs(diff(pgf, x(i)),x,ones(1,q))); 
end
z=expect_pgf;
end

% variance from pgf: second factorial moment=variance-expectation =second derivative of pgf
% valued at 1 - first derivative of pgf valued at 1 
function [z]=Var_pgf(pgf,x,q,qi)   
Var_pgf=zeros(q^2,qi);
for i=1:q
for j=1:q
store(i,j,:)=(j~=i)*((subs(diff(diff(pgf, x(j)),x(i)),x,ones(1,q)))-(subs(diff(pgf, x(j)),x,ones(1,q))).*(subs(diff(pgf, x(i)),x,ones(1,q)))); %COV+EiEj-EiEj=COV(i,j) pour i diff j
end
end
for i=1:q
store2(i,:)=(subs(diff(diff(pgf, x(i)),x(i)),x,ones(1,q)))-(subs(diff(pgf, x(i)),x,ones(1,q))).^2 + (subs(diff(pgf, x(i)),x,ones(1,q))); %Var+EiEi-Ei - EiEi+Ei = Var?
end
for i=1:qi
Covar= store(:,:,i)+diag(store2(:,i));
Var_pgf(:,i)=Covar(:);
end
z=Var_pgf;
end

% distribution of univariate r.v. from pgf, by successive derivations and
% valuations at 0 
function [z]=dist_pgf(pgf,y,maxn)  
Dist=zeros(2,maxn);
Der_n=pgf;
for n=0:maxn
Dist(1,n+1)=n;
Dist(2,n+1)=(1/factorial(n))*subs(Der_n,0);
Der_n=diff(Der_n);
end
z=Dist;
end
