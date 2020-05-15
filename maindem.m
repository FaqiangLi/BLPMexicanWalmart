% This program computes the random coefficeints discrete choice model
% described in thhe "A Research Assistant's Guide to Discrete Choice Models
% of Demand," NBER technical paper #221, and "Measuring Market Power in 
% the Ready-to-Eat Cereal Industry," NBER WP #6387.

% Written by Aviv Nevo, May 1998.  
% Updated to matlab 6 and 7 by Eric Rasmusen, Dec. 31, 2005 and 
% Bronwyn Hall, April 2005 

%  this program does not include the observable consumer demographics.
%  Use the file main.m, rather than this one, to include observable demog.

%The inputs to this code are the data matrices ps2.mat and iv.mat
%The output from this program will be results.txt and  mydiary.txt.  

disp('                                                                   ')
disp('*******************************************************************')
disp('                                                                   ')

delete mydiary.txt % erase any old mydiary.txt file.
diary mydiary.txt

date = datestr(clock)

global invA ns x1 x2 s_jt IV vfull dfull theta1 theti thetj cdid cdindex
%global nmkt nbrn n_inst

% load data. see description in data_readme.pdf
% ps2.mat contains the matrices/vectors id, v, demogr, x1, x2,s_jt, id_demo

load ps2.mat
load iv.mat

disp('                                                                  ')
disp('******************************************************************')
disp('                                                                  ')
ns = 20;       % number of simulated "individuals" per market %
nmkt = 94;     % number of markets = (47 cities)*(2 quarters)  %
nbrn = 24;     % number of brands per market. 
               % We have 24*94=2256 observations %
n_inst = 20 ;  % Number of instruments for price (other time periods)%

disp('                                                                  ')
disp('******************************************************************')
disp('                                                                  ')

% Creates a matrix of instruments from the last 20 cols of iv.mat and the 
% 24 brand dummies contained in x1

IV = [iv(:,2:(n_inst+1)) x1(:,2:(nbrn+1))];

disp(['The dimensions of object IV are ' num2str(size(IV)) ])

clear iv

disp('                                                                   ')
disp('*******************************************************************')
disp('                                                                   ')

% Create a vector identifying the market for each observation.
% cdid = [1...1,2...2,3...3,...94...94]'

cdid = kron([1:nmkt]',ones(nbrn,1));
disp(['Object cdid dimensions: ' num2str(size(cdid)) ])

% Create a vector identifying the last brand observation for each market
% Each block of 24 observations is for a single market
% cdindex = [24 48 72....2256]'

cdindex = [nbrn:nbrn:nbrn*nmkt]';
disp(['The dimensions of object cdindex are     ' num2str(size(cdindex)) ])
disp('                                                                   ')
disp('*******************************************************************')
disp('                                                                   ')

% starting values. These are the converged values from Nevo's JEMS paper
% Coeff that will not be maxed over are fixed at zero.   %
% rows correspond to constant, price, sugar, mushy
% columns correspond to sigma, income, income^2, age, child
% There are 13 parameters to be estimated, and 7 set to zero by assumption
% These are 4 characteristic parameters and 9 interaction parameters.
% Some interactions are set to zero (e.g. Sugar and Income squared)

%theta2w=    [0.5      2.0           0      1.5           0;
 %            3.0     10.1        -0.2        0         6.5 ;
 %           0.02   -0.25           0      0.5           0;
 %            0.50    2.3            0     -2.0           0];
         
% starting values for the nonlinear parameters in theta2  %
% rows correspond to constant, price, sugar, mushy
% columns correspond to sigma, income, income^2, age, child
% There are 4 parameters to be estimated one for each characteristic

theta2w=    [0.5  ;
             3.0  ;
             0.02 ;
             0.50 ];

disp('                                                                   ')
disp(['The dimensions of object theta2w are     ' num2str(size(theta2w)) ])
disp('                                                                   ')

% create a vector of the non-zero elements in the above matrix, and the %
% corresponding row and column indices. this facilitates passing values %
% to the functions below. %

[theti, thetj, theta2]=find(theta2w) 

horz=['    mean       sigma    '];
vert=['constant  ';
      'price     ';
      'sugar     ';
      'mushy     '];

disp('                                                                   ')
disp('*******************************************************************')
disp('                                                                   ')

	% create weight matrix

invA = inv([IV'*IV]);

disp(['The dimensions of object invA  are ' num2str(size(invA)) ])
disp('                                                                   ')
disp('*******************************************************************')
disp('                                                                   ')

% Simple Logit results - each brand share is normalized by the share of 
% the outside good in that market.
% save the mean utilities for each brand/mkt to use 
% as initial values for the search below.

temp = cumsum(s_jt);
sum1 = temp(cdindex,:);
sum1(2:size(sum1,1),:) = diff(sum1);
outshr = 1.0 - sum1(cdid,:);

disp(['Object temp dimensions: ' num2str(size(temp)) ])
disp(['Object sum1 dimensions: ' num2str(size(sum1)) ])
disp(['Object outshr dimensions: ' num2str(size(outshr)) ])

y = log(s_jt) - log(outshr);
mid = x1'*IV*invA*IV';
t = inv(mid*x1)*mid*y;  %  price and brand dummies only
mvalold = x1*t;
oldt2 = zeros(size(theta2));
mvalold = exp(mvalold);

disp(['Object y dimensions: ' num2str(size(y)) ])
disp(['Object mid dimensions: ' num2str(size(mid)) ])
disp(['Object t dimensions: ' num2str(size(t)) ])
disp([' Simple logit estimates'])
t
disp(['Object mvalold dimensions: ' num2str(size(mvalold)) ])
disp(['Object oldt2 dimensions: ' num2str(size(oldt2)) ])

% the next command creates a new file, mvalold.mat, with mvaold 
% and oldt2 in it, and then clears the logit information from memory. 
% disp(['   mvalold from simple logit'])
% mvaold=mvalold(1:30)

save mvalold mvalold oldt2
clear mid y outshr t oldt2 mvalold temp sum1

%  set up the consumer data for each brand/mkt.  Replicate it for each
% brand in the market

vfull = v(cdid,:);
% dfull = demogr(cdid,:);

disp(['Object vfull dimensions: ' num2str(size(vfull)) ])
disp(['Object dfull dimensions: ' num2str(size(dfull)) ])
disp('                                                                   ')
disp('*******************************************************************')
disp('                                                                   ')

%   GMM estimation

mytolfun=0.000001 ;
mytolx = 0.000001 ;
mytoliter = 5000 ;
mytolfeval= 5000 ;

%{
 method(1) Computes the estimates using a Quasi-Newton method
 with an *analytic* gradient. Fminunc is from "optimization toolbox" 
%}

options=optimset('GradObj','on','HessUpdate','bfgs','LargeScale','off',...
    'MaxIter',mytoliter,'Display','iter','TolFun',mytolfun,...
    'TolX',mytolx,'DerivativeCheck','off','MaxFunEvals',mytolfeval)
tic
%[theta2,fval,exitflag,output,grad,hess] = fminunc('gmmobjg',theta2,options)
[theta2,fval,exitflag,output,grad,hess] = fminunc('gmmobjgdem',theta2,options)
comp_t = toc/60;

%{
 method(2) computes the estimates using a simplex search method.

options=optimset('MaxIter',mytoliter,'Display','iter','TolFun',mytolfun,...
    'TolX',mytolx,'MaxFunEvals',mytolfeval)
tic  
[theta2,fval,exitflag,output]  = fminsearch('gmmobjnmdem',theta2, options)
comp_t = toc/60 ;
%}

disp('                                                                  ')
disp('*******************************************************************')
disp('                                                                   ')

	% computing the s.e.
    %var_cov is a separate matlab file

vcov = var_covdem(theta2);
se = sqrt(diag(vcov));
se

disp(['Object vcov dimensions: ' num2str(size(vcov)) ])
disp(['Object se dimensions: ' num2str(size(se)) ])

theta2w = full(sparse(theti,thetj,theta2));
theta2w
t = size(se,1) - size(theta2,1);
se2w = full(sparse(theti,thetj,se(t+1:size(se,1))));
se2w

disp(['Object theta2w dimensions:     ' num2str(size(theta2w)) ])
disp(['Object t dimensions:     ' num2str(size(t)) ])
disp(['Object se2w dimensions:     ' num2str(size(se2w)) ])
disp('                                                                  ')
disp('******************************************************************')
disp('                                                                  ')

	% the MD estimates

omega = inv(vcov(2:25,2:25));
xmd = [x2(1:24,1) x2(1:24,3:4)];
ymd = theta1(2:25);
theta1

disp(['Object omega dimensions: ' num2str(size(omega)) ])
disp(['Object xmd dimensions: ' num2str(size(xmd)) ])
disp(['Object ymd dimensions: ' num2str(size(ymd)) ])

beta = inv(xmd'*omega*xmd)*xmd'*omega*ymd;
beta
resmd = ymd - xmd*beta;
semd = sqrt(diag(inv(xmd'*omega*xmd)));
semd
mcoef = [beta(1); theta1(1); beta(2:3)];
mcoef
semcoef = [semd(1); se(1); semd(2:3)];
semcoef

disp(['Object omega dimensions: ' num2str(size(omega)) ])
disp(['Object beta dimensions: ' num2str(size(beta)) ])
disp(['Object resmd dimensions: ' num2str(size(resmd)) ])
disp(['Object semd dimensions: ' num2str(size(semd)) ])
disp(['Object mcoef dimensions: ' num2str(size(mcoef)) ])
disp(['Objectsemcoef dimensions: ' num2str(size(semcoef)) ])

Rsq = 1-((resmd-mean(resmd))'*(resmd-mean(resmd)))/...
                     ((ymd-mean(ymd))'*(ymd-mean(ymd)));
Rsq_G = 1-(resmd'*omega*resmd)/((ymd-mean(ymd))'*omega*(ymd-mean(ymd)));
Chisq = size(id,1)*resmd'*omega*resmd;

disp('                                                                   ')
disp('*******************************************************************')
disp('                                                                   ')

% diary results.txt  %  this only contains the parameter est and se.  
disp(horz)
disp('  ')
for i=1:size(theta2w,1)
     disp(vert(i,:))
     disp([mcoef(i) theta2w(i,:)])
     disp([semcoef(i) se2w(i)])
end

disp('                                                                   ')
disp('*******************************************************************')
disp('                                                                   ')
disp(['GMM objective:  ' num2str(fval)])
disp(['MD R-squared:  ' num2str(Rsq)])
disp(['MD weighted R-squared:  ' num2str(Rsq_G)])
disp(['run time (minutes):  ' num2str(comp_t)])
disp('                                                                   ')
disp('*******************************************************************')
disp('                                                                   ')

diary off
delete gmmresid.mat 
delete mvalold.mat
