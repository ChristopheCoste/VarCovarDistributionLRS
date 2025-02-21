clear all; clc; close all;

%% inputs
epsilon=1e-14 % tolerance for approximations (=epsilon1=epsilon2)
% expectation of survival
S=[0.9030  0 0  0 0 0 ;
   0.0038 0.96070 0 0 0 0;
   0 0.01225 0.96545 0 0 0 ;
   0 0 0.01735 0.97595 0 0;
   0 0 0  0.01205 0.96335 0;
   0 0 0 0 0.01835 0.9903]
s=size(S,1) % number of classes

F=zeros(s,s);%expectation of reproduction
F(1,:)=[0 0 0.299 0.77415 1.9573 6.0251]

%% probability distribution of (time-step) reproduction
epsilon2=epsilon;
% reproduction is poisson, stage 6 has highest parameter
% we truncate it according to epsilon2
f6=F(1,6);
Pf6=[0];
i=0;
while  max(1-sum(Pf6))>epsilon2
Pf6(i+1,:)=(f6^i)*(exp(-f6))*1/factorial(i);
i=i+1;
end
alpha=length(Pf6)-1 %maximum number of offspring produced per time-step considered

Fcal=zeros(alpha+1,s); % we truncate the reproduction of all stages according to N
for j=1:s
for i=0:alpha
Fcal(i+1,j)=(F(1,j)^i)*(exp(-F(1,j)))*1/factorial(i);
end
end

Fcal

%% maximum age
epsilon1=epsilon;
% the model has no maximum age, we compute an approximation for it so that
% the probability to survive past this age is below epsilon
omega=0;
smax=1; 
while smax>epsilon
omega=omega+1;
smax=max(ones(1,s)*(S^omega));
end
omega 

%% probability distribution of the number of newborn

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
Rcal;



%% asymptote
probgrowth=F(1,6)/(F(1,6)-log(S(6,6)))
Rcal(1001,:)./Rcal(1000,:)

%% figures 
figure
t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');

nexttile;
 kmax=20; %max number for figure
semilogy(0:kmax,(Rcal(1:kmax+1,1)),'LineWidth',4)
hold on
semilogy(0:kmax,(Rcal(1:kmax+1,2:end)),'LineWidth',2)
xlabel('lifetime reproductive success','Interpreter','latex')
ylabel('probability','Interpreter','latex')
legend('LRS','event. num. of offs. prod. by stage 2','event. num. of offs. prod. by stage 3','event. num. of offs. prod. by stage 4','event. num. of offs. prod. by stage 5','event. num. of offs. prod. by stage 6','Interpreter','latex')

kmax=10000; %max number for figure
kmax=2000; %max number for figure

nexttile;
semilogy(0:kmax,(Rcal(1:kmax+1,1)),'LineWidth',4)
hold on
semilogy(0:kmax,(Rcal(1:kmax+1,2:end)),'LineWidth',2)
hold on
semilogy(0:kmax,probgrowth.^(0:kmax),'LineWidth',4)
xlabel('lifetime reproductive success','Interpreter','latex')
ylabel('probability','Interpreter','latex')
legend('LRS','event. num. of offs. prod. by stage 2','event. num. of offs. prod. by stage 3','event. num. of offs. prod. by stage 4','event. num. of offs. prod. by stage 5','event. num. of offs. prod. by stage 6','asymptote $\rho^{k}$','Interpreter','latex')

figure
kmax=1000
loglog(4:kmax,Rcal(4:kmax,1)./Rcal(3:kmax-1,1),'LineWidth',4)
hold on
loglog(4:kmax,Rcal(4:kmax,2:end)./Rcal(3:kmax-1,2:end),'LineWidth',2)
hold on
loglog(4:kmax,probgrowth*ones(kmax-3,1),'LineWidth',4)

ylabel('Probability Progression Ratio $\frac{P_{i+1}}{P_{i}}$','Interpreter','latex')
xlabel('lifetime reproductive success','Interpreter','latex')
legend('LRS','event. num. of offs. prod. by stage 2','event. num. of offs. prod. by stage 3','event. num. of offs. prod. by stage 4','event. num. of offs. prod. by stage 5','event. num. of offs. prod. by stage 6','asymptote $\rho$','Interpreter','latex')
ylim([.9 1.3])

%% probability to reach stage 3 from stage 1
S3=S(1:3,1:3);
S3(3,3)=0;
N3=(eye(3)-S3)^-1;
N3(3,1)
Rcal(1:10,1)./Rcal(1:10,3)