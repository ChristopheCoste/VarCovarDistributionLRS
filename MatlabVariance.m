clear all; clc;close all;

% Provides the expectation and variance-covariance of LRS
% Can be used for any strutured population projeciton model
% We use the model of the illustration.

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

%% projection of abundance vector over time
n0=[1;1]
Vn1=(VS+VF)*n0
reshape(Vn1,s,s)


%% expectation of LRS
R=F*(eye(s)-S)^(-1) %eq 8, next generation matrix

rT=ones(1,s)*F*(eye(s)-S)^(-1) % expectation of Total LRS



%% Variance-Covariance of LRS
N=(eye(s)-S)^(-1); % fundamental matrix: expected time-spent in stages

VR= (VF+kron(R,R)*VS)*N %eq9

reshape(VR(:,1),s,s) %var-covariance of LRS for a type 2 (born medium) individual


%% variance of total LRS

VrT=(kron(rT,rT)*VS+ones(1,s^2)*VF)*(eye(s)-S)^(-1) %eq 11

