clear all; clc;close all;

%% inputs

% number of states
s=2 

%expectation of survival
S=[.3 0; .2 .7] % expectation of survival

%expectation of reproduction
F=[0.5 0.5; 0 0.5] 

%variance of reproduction
VF=[.25 .25; 0 -.25 ; 0 -.25 ;0 0.25]

%% intermediary outputs
% vector of expectation of reproduction 
fT=ones(1,s)*F

% variance survival
J=zeros(s^2,s);
for i=1:s
ei=zeros(s,1);ei(i)=1;
J=J+kron(ei,ei)*ei';
end

VS=J*S-kron(S,S)*J    %(eq.4)



%% expectation of LRS
R=F*(eye(s)-S)^(-1) %eq 8, next generation matrix

rT=ones(1,s)*F*(eye(s)-S)^(-1) % expectation of Total LRS


%% covariance between survival and fert
C=zeros(s^2,s);

C(1,1)=-0.15;
C(2,1)=-0.1;
C

%%
K=zeros(s^2,s^2)
for i=1:s
    for j=1:s
ei=zeros(s,1); ej=zeros(s,1);   ei(i)=1; ej(j)=1;
K=K+kron(ei,ej)*kron(ej,ei)';
    end
end
K

%% variance of LRS

N=(eye(s)-S)^(-1); % fundamental matrix: expected time-spent in stages


VR=(VF+kron(R,R)*VS)*N + (kron(eye(s),R) + kron(R,eye(s))*K)*C*N 



%% alternative method

%inputs
S=[0 .3 0; 0 .3 0;0 .4 .7] 
s=size(S,1) 
F=[0.5 0 0.25;0.5 0 0.25; 0 0 0.5] 

% variance of F and S
VF=BernV(F)
VS=BernV(S) 

%
N2=(eye(s)-S)^(-1); 
R2=F*N2; 
VR2= (kron(R2,R2)*VS+VF)*N2 

%%
P1=[0.5 0; 0.5 0;0 1]
P2=[1 1 0; 0 0 1]

R_2=P2*R2*P1

VR_2=kron(P2,P2)*(kron(R2,R2)*BernV(P1)+VR2*P1)

%% functions
function [ V] = BernV(A)  % eq 35, Function providing the variance of a Bernoulli scheme
s=size(A,2);
a=size(A,1);
V=zeros(a*a,s);
for i=1:s
vecc=A(:,i);
V(:,i)=vec(diag(vecc)-vecc*vecc');
end
end

function [y] = vec(x)
    y=x(:);
end