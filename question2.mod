// Representantive household with EZ preferences
//
// Constanza Vergara

//----------------------------------------------------------------
// 0. Housekeeping
//----------------------------------------------------------------

close all

//----------------------------------------------------------------
// 1. Endogenous variables (17=3+5+4+1+2+1+1)
//----------------------------------------------------------------

var 

// Utility (3)
v v1 dv1

// Allocation (5)
y c i l k

// Prices (4)
r w lambda m

// Preferences (1)
n

// Taxes (2)
taul tauk

// Government Expenditure (1)
g

// Productivity (1)
z;

//----------------------------------------------------------------
// 2. Exogenous variables (4)
//----------------------------------------------------------------

varexo 

ez en el ek;

//----------------------------------------------------------------
// 3. Parameters
//----------------------------------------------------------------

parameters 

// Adjustment to the utility function
cte

// Utility Function
ggamma ttheta bbeta

// Preferences
rrho_n ssigma_n mmean_n

// Taxes
rrho_l rrho_k ssigma_l ssigma_k mmean_tl mmean_tk

//Productivity
rrho_z ttheta_z1 ttheta_z2 aalpha ddelta ppsi;


//----------------------------------------------------------------
// 4. Calibration
//----------------------------------------------------------------

// Adjustment to the utility function
 cte=3;

//Utility Function
 ggamma=10;
 ttheta=-18;
 bbeta=.99;

//Preferences
 rrho_n=.95;
 ssigma_n=.002;
 mmean_n=1;

// Taxes
 rrho_l=.95;
 ssigma_l=0.005;
 rrho_k=.9;
 ssigma_k=.008;
 mmean_tl=.25;
 mmean_tk=.35;

// Productivity
 rrho_z=.95;
 ttheta_z1=.007;
 ttheta_z2=-.001; 
 aalpha=.33;
 ddelta=.03;
 ppsi=.05;


//----------------------------------------------------------------
// 5. Computation of Steady State
//----------------------------------------------------------------
n_ss=mmean_n;
taul_ss=mmean_tl;
tauk_ss=mmean_tk;
z_ss=0;

r_ss=(1/bbeta-(1-ddelta))/(1-tauk_ss);
w_ss=(1-aalpha)*(aalpha/r_ss)^(aalpha/(1-aalpha));
l_ss=sqrt((1-taul_ss)*w_ss/(n_ss*((aalpha/r_ss)^(aalpha/(1-aalpha))*(1-(1-aalpha)*taul_ss-aalpha*tauk_ss)-ddelta*(aalpha/r_ss)^(1/(1-aalpha)))));

k_ss=l_ss*(aalpha/r_ss)^(1/(1-aalpha));
i_ss=ddelta*k_ss;
g_ss=taul_ss*w_ss*l_ss+tauk_ss*r_ss*k_ss;
y_ss=k_ss^aalpha*l_ss^(1-aalpha);
c_ss=y_ss-i_ss-g_ss;

v_ss=(log(c_ss)-n_ss*l_ss^2/2+cte);
v1_ss=v_ss^(1-ggamma);
dv1_ss=v_ss^-ggamma;
lambda_ss=(1-bbeta)/c_ss;
m_ss=lambda_ss;


//----------------------------------------------------------------
// 6. Model
//----------------------------------------------------------------

model; 

 // 1. Value function
  v = ((1-bbeta)*(log(c)-n*l^2/2+cte)^((1-ggamma)/ttheta)+bbeta*v1^(1/ttheta))^(ttheta/(1-ggamma));

 // 2. Value Function to the power of -9
 v1=v(+1)^(1-ggamma);

 // 3. Derivative of the value function to the power of -9
 dv1=v(+1)^(-ggamma);

 // 4. Preferences
 (n-mmean_n)=rrho_n*(n(-1)-mmean_n)+ssigma_n*en;

 // 5. Labor Tax
 (taul-mmean_tl)=rrho_l*(taul(-1)-mmean_tl)+ssigma_l*el;

 // 6. Capital Tax
 (tauk-mmean_tk)=rrho_k*(tauk(-1)-mmean_tk)+ssigma_k*ek;

 // 7. Rental Rate capital
 r=exp(z)*aalpha*(l/k(-1))^(1-aalpha);

 // 8. Wages
 w=exp(z)*(1-aalpha)*(k(-1)/l)^(aalpha);

 // 9. Law of motion of capital
 k=(1-ddelta)*k(-1)+(1-ppsi*(i/i(-1)-1)^2)*i;

 // 10. Government Expenditure
 g=taul*w*l+tauk*r*k(-1);

 // 11. Production Function
 y=exp(z)*k(-1)^(aalpha)*l^(1-aalpha);

 // 12. Resources Constraint of the Economy
 c=y-i-g;

 // 13. Productivity Process
 z=rrho_z*z(-1)+ttheta_z1*ez+ttheta_z2*ez(-1);

 // 14. Static Condition leisure-consumption
 (-n*l+(1-taul)*w/c)=0;

 // 15. FOC with respect to consumption
lambda=v^((ttheta-1+ggamma)/ttheta)*(1-bbeta)/c*(log(c)-n*l^2/2+cte)^((1-ggamma-ttheta)/ttheta);

// 16. FOC with respect to k prime
m=bbeta*v^((ttheta-1+ggamma)/ttheta)*v1^(1/ttheta-1)*dv1*(lambda(+1)*r(+1)*(1-tauk(+1))+(1-ddelta)*m(+1));

// 17. FOC with respect to investment
lambda-m*(1-ppsi*(i/i(-1)-1)^2-2*ppsi*(i/i(-1)-1)*(i/i(-1)))=
bbeta*v^((ttheta-1+ggamma)/ttheta)*v1^(1/ttheta-1)*dv1*m(+1)*2*ppsi*(i(+1)/i-1)*(i(1)/i)^2;

end;

initval;
    v=v_ss;
    v1=v1_ss;
    y=y_ss;
    c=c_ss;
    i=i_ss;
    l=l_ss;
    k=k_ss;
    r=r_ss;
    w=w_ss;
    n=n_ss; 
    taul=taul_ss;
    tauk=tauk_ss;
    g=g_ss;
    z=z_ss;
    m=m_ss;
    lambda=lambda_ss;
    dv1=dv1_ss;

    ez=0;
    en=0;
    el=0;
    ek=0;

end;

shocks;
  var ez=1;
  var en=1 ;
  var el=1; 
  var ek=1 ;
end;

steady;

stoch_simul(hp_filter = 1600, irf = 100, order = 3);

