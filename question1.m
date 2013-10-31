%% Homework 3 Econ 714
%
% Chebyshev
% Constanza Vergara

%% 0. Housekeeping

clear all
close all
clc

tic

%%  1. Calibration

aalpha = 1/3;           % Elasticity of output w.r.t. capital
bbeta  = 0.95;          % Discount factor
ddelta = .09;           % Depreciation
laborSteadyState=1/3;   % Labor in Steady State

%Grid Points
nGridCapital = 1728;

% Productivity values
vProductivity = [exp(-0.04);exp(0);exp(0.04)]';

% Transition matrix
mTransition   = [0.9727, 0.0273, 0.0000;
                 0.0041, 0.9806, 0.0153;
                 0.0000, 0.0082, 0.9836];
             
%% 2. Steady State

capitalSteadyState = ((1/bbeta-1+ddelta)/(aalpha*laborSteadyState^(1-aalpha)))^(1/(aalpha-1));
outputSteadyState = capitalSteadyState^aalpha*laborSteadyState^(1-aalpha);
consumptionSteadyState = outputSteadyState-capitalSteadyState*ddelta;

psi=outputSteadyState*(1-aalpha)/((outputSteadyState-ddelta*capitalSteadyState)*laborSteadyState^2);
nGridProductivity = length(vProductivity);

vGridCapital = linspace(0.5*capitalSteadyState,1.5*capitalSteadyState,nGridCapital);

%% 3. Chebyshev Matrices


order=8;
theta0=zeros(order*nGridProductivity,1);

a=0.5*capitalSteadyState;
b=1.5*capitalSteadyState;

x=col_points(a,b,order);

myfun=@(theta,x,a,b) systEquation(theta,x,vProductivity,mTransition,aalpha,psi,bbeta,ddelta,a,b);

for i=0:2
    theta0(order*i+1)=laborSteadyState;
    theta0(order*i+2)=-0.04;
    theta0(order*i+3)=.01;
end   

options=optimset('Display','off');
[theta,fval] = fsolve(myfun,theta0,options,x,a,b);

theta=real(theta);

mLaborFunction=zeros(nGridCapital,nGridProductivity);

for nProductivity=1:nGridProductivity
    mLaborFunction(:,nProductivity)=cheby_approx(vGridCapital,theta((nProductivity-1)*order+1:...
        nProductivity*order)',a,b);
end

mPolicyFunction=zeros(nGridCapital,nGridProductivity);

for nProductivity=1:nGridProductivity
    mPolicyFunction(:,nProductivity)=vGridCapital'*(1-ddelta)+vProductivity(nProductivity).*...
        vGridCapital'.^aalpha.*mLaborFunction(:,nProductivity).^(1-aalpha)-...
        (consumption(vProductivity(nProductivity),...
        vGridCapital,aalpha,psi,mLaborFunction(:,nProductivity)'))';
end

residual=zeros(1,nGridProductivity);

for nProductivity=1:nGridProductivity
    residual(nProductivity)=mean(abs(euler2(vProductivity(nProductivity),vGridCapital,vProductivity,...
        mTransition,aalpha,psi,bbeta,ddelta,theta,a,b)));
end


%% 6. Plotting results

figure(1)

subplot(2,1,1)
plot(vGridCapital,mPolicyFunction)
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
title('Policy Function')

subplot(2,1,2)
plot(vGridCapital,mLaborFunction)
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
title('Labor Function')

toc
