%% inputs
%tolerance 
epsilon1=0.000000000001 %proportion of individuals allowed to survive past the maximum age
s=2 % number of states
S=[.3 0; .2 .7] % expectation of survival
F=[0.5 0.5; 0 0.5] %expectation of reproduction

%numeric matrix storing (in F(i,j) the probability distribution for an individual in state j to produce i-1 offspring in a given time-step
Fcal=[.5 0; .5 1]


alpha=size(Fcal,1)-1  %maximum number of offspring that can be produced by an individual in a time step


%% maximum age
% the model has no maximum age, we compute an approximation for it so that
% the probability to survive past this age is below epsilon
omega=0;
smax=1; 
while smax>epsilon1
omega=omega+1;
smax=max(ones(1,s)*(S^omega));
end
omega 

%% generation of probability distribution of LRS by iteration
maxr=alpha*omega  %maximum coeff of polynomial of r = maximum LRS allowed by the approximation+ 1

Rcal=zeros(maxr+1,s);
Rcal(1:alpha+1,:)=Fcal; %initialisation (first line of equation 18)

for age=omega-1:-1:1 %iteration of equation 18
Qcal=Rcal*S ;
Qcal(1,:)=Qcal(1,:)+1-sum(S);
for j=1:s   % convolution of equation 18
vect= flip(conv(flip(Qcal(:,j)),flip(Fcal(:,j)))); 
Rcal(:,j)=vect(1:maxr+1);
end
end

%% figure
kmax=15; %max number for figure
 figure
 t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
nexttile;
bar(0:kmax,Rcal(1:kmax+1,1),'b')
xlabel('lifetime reproductive success','Interpreter','latex')
ylabel('probability','Interpreter','latex')
legend('total LRS for a type \emph{small}','Interpreter','latex')
ylim([0 0.4])

nexttile;
bar(0:kmax,Rcal(1:kmax+1,2),'r')
xlabel('lifetime reproductive success','Interpreter','latex')
ylabel('probability','Interpreter','latex')
legend('total LRS for a type \emph{small}','Interpreter','latex')
ylim([0 0.4])
