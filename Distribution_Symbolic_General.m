clear all; clc;close all;

%% input1: expectation of survival
%there are 3 states that are types --> S is diagonal
S=[.6 0 0; 0 .5 0 ;0 0 .7] % expectation of survival

%% number of states and dummy variables (output)
s=size(S,1) % number of states

x = sym('x_%d',[1 s]);
y = sym('y');

%% input2: pgf of reproduction

pgfF1= 0.5*x(1) + 0.5 % type 1 individual produce one type 1 individual with probability 0.5
pgfF2= exp(x(1)-1)*exp(x(2)-1) % type 2 individuals produce type 1 and 2 individuals according to 2 independent Poisson of parameter 1
pgfF3= 0.5*x(1) + 0.2*x(3)+0.3 % type 3 individuals produce 0 offspring with prob .3, a type 1 with prob .5 and a type 3 with probability .2
pgfF  =[pgfF1, pgfF2, pgfF3]

%% intermedidary output: expectation of reproduction
F=expect_pgf(pgfF,x,s,s) % for sanity check

%% intermedidary output: pgf of survival 
pgfS  = pgfBern([1:s],S,x); %for sanity check
vpa(pgfS,2)

%% OUTPUT1: pgf of R
dd=diag(pgfF);
eye(s)-S*dd
invinv=(eye(s)-S*dd)^(-1)
pgfR=(pgfF-ones(1,s)*S*dd)*invinv %
vpa(pgfR,2)

% pgf of rT
pgff=subs(pgfF,x,y*ones(1,s)); %pgf of univariate (i.e. unstructured) reproduction
dd=diag(pgff);
pgfr=(pgff-ones(1,s)*S*dd)*(eye(s)-S*dd)^(-1); %eq 14b
vpa(pgfr,2)

%% OUTPUT2: moments of LRS 
expect_pgf(pgfR,x,s,s)
Var_pgf(pgfR,x,s,s)
sum(Var_pgf(pgfR,x,s,s))

expect_pgf(pgfr,y,1,s)
Var_pgf(pgfr,y,1,s)


%% MAIN OUTPUT: joint probability distribution of LRS , obtained by succesive derivations of the pgfs

typepar=2; %Choose the type of the parent
LRS=[100;50;0]; % chose the LRS (a vector of s integers) of which one wants to know the probability


Prob=1;
Der=pgfR(typepar);
for i=1:s
Prob=Prob*(1/factorial(LRS(i)));
for ni=1:LRS(i)
Der=diff(Der,x(i));
end
end
Der;

Prob=Prob*subs(Der,{x(1),x(2),x(3)},{0,0,0});
%vpa(Prob,2)

fprintf('The probability that an individual in state %g\n eventually produces this vector of offspring :',typepar)
LRS
fprintf(' is %g\n',Prob)
fprintf('If state %g is a type, the eventual number of offspring is the LRS',typepar)


%% functions

%function providing the pgf of a Bernoulli scheme (from equation 36)    
function [z]=pgfBern(i,A,x)                     
z = 1+(x-1)*A(:,i);
end

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
