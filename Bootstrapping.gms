*****Scenario Generation*****

***Loading Data file.***

SETS
         Date            'Dates'
         Asset           'Assets'
         Weeks           'Weeks for bootstrapping' /w1*w4/
         Scenarios       'Scenarios for bootstrapping' /s1*s250/
;
ALIAS(Date,t);
ALIAS(Asset,i);
ALIAS(Weeks,w);
ALIAS(Scenarios,s,l);


PARAMETER
         Return(Date,Asset)              'Return of index i at date t.'
         ScenarioReturnW(i, s, w)        'Return of index i, at scenario s, in week w'
         ScenarioReturn(i,s)             'Return of index i, at scenario s'
         oneAsset(t)
;

SCALARS
         RandomNumber    'A random number for bootstrapping';

*Loading date,asset and ETF returns data from "out-ETFRet.gdx"
$GDXIN out-ETFRet
$LOAD  Date,Asset,Return
$GDXIN



*Subset of 31 ETFs as a result from Minimum Spanning Tree method
SET
    SubAsset(Asset) /BRF,DBA,DBB,EIDO,EIRL,ENZL,FXB,FXE,FXF,IAU,IGN,IGV,IHI,IPFF,ITB,IYZ,MLPA,PALL,PGF,PGX,SOCL,TAN,TUR,UGA,URA,USDU,VNM,XBI,XHS,XMLV,XPH/
;
display SubAsset;

oneAsset(t) = Return(t , 'ANGL');
display t, i, Return, oneAsset;



PARAMETER
   scenario_ids(Weeks, Scenarios) All returns dataframe;





LOOP(s,
    loop(w,
        RandomNumber = uniformint(1,97);
        ScenarioReturnW(i, s, w) = sum(t$(ord(t)=RandomNumber), Return(t,i));
     );
     ScenarioReturn(i,s) =  sum(w, ScenarioReturnW(i, s, w));
);
EXECUTE_UNLOAD 'FixedScenarios.gdx', ScenarioReturn;
EXECUTE 'gdxxrw.exe FixedScenarios.gdx O=FixedScenarios.xls par=ScenarioReturn rng=sheet1!a1' ;

display ScenarioReturn;

***********************QUESTIN 2**########




parameter WorstScen(i);
WorstScen(i) =  smin(l, ScenarioReturn ( i, l ));
display WorstScen;

parameter SumScenRet(l), WC, BC, AC;
SumScenRet(l) = sum(i, ScenarioReturn ( i, l ))/card(i);
WC=smin(l, SumScenRet(l));
BC=smax(l, SumScenRet(l));
AC=sum(l,SumScenRet(l))/card(l);
display SumScenRet, WC, BC, AC, ScenarioReturn;


SCALARS
        Budget        'Nominal investment budget'
        alpha         'Confidence level'
        MU_Target     'Target portfolio return'
        MIN_MU        'Minimum return in universe'
        MAX_MU        'Maximum return in universe'
;

Budget = 100.0;
alpha  = 0.99;

PARAMETERS
        pr(l)       'Scenario probability'
        P(i,l)      'Final values'
        EP(i)       'Expected final values'
;

pr(l) = 1.0 / CARD(l);

P(i,l) = 1 + ScenarioReturn ( i, l );

EP(i) = SUM(l, pr(l) * P(i,l));

MIN_MU = SMIN(i, EP(i));
MAX_MU = SMAX(i, EP(i));
display MIN_MU, MAX_MU, pr;

MU_TARGET = MIN_MU;
MU_TARGET = MAX_MU;
display MU_TARGET;




*From here

scalar HighestLoss, StepSize;
SCALAR
     lambda    'Risk attitude'
;
HighestLoss = Budget*(smax((i,l), P(i,l) )- smin((i,l), P(i,l)));
display P, HighestLoss,Budget;


POSITIVE VARIABLES
        x(i)            'Holdings of assets in monetary units (not proportions)'
        VaRDev(l)       'Measures of the deviations from the VaR'
;

VARIABLES
       VaR             'Value-at-Risk'
       CVaR            'Objective function value'
       Losses(l)       'Measures of the losses'
       CVARLambda
;

Binary variable
      Z(l)             '1 if VarDevCon(l)>0 otherwise 0'

EQUATIONS
        BudgetCon        'Equation defining the budget contraint'
        ReturnCon        'Equation defining the portfolio return constraint'
        LossDef(l)       'Equations defining the losses'
        VaRDevCon(l)     'Equations defining the VaR deviation constraints'
        ObjDefCVaR       'Objective function definition for CVaR minimization'
        SetZ(l)
        ObjDefVaR
        ObjDefLambda
;

BudgetCon ..         SUM(i, x(i)) =E= Budget;

ReturnCon ..         SUM(i, EP(i) * x(i)) =G= MU_TARGET * Budget;

LossDef(l)..         Losses(l) =E= (Budget - SUM(i, P(i,l) * x(i)));
*LossDef(l)..         Losses(l) =E= - SUM(i, P(i,l) * x(i));

VaRDevCon(l) ..      VaRDev(l) =G= Losses(l) - VaR;

ObjDefCVaR ..        CVaR =E= VaR + SUM(l, pr(l) * VaRDev(l)) / (1 - alpha);

ObjDefLambda    ..   CVARLambda     =E= (1-lambda) * SUM(i, EP(i) * x(i)) - lambda * CVaR;

SetZ(l)  ..          HighestLoss*z(l) =G= VaRDev(l);

ObjDefVaR ..         1 - sum(l, pr(l)*z(l)) =G= alpha;

MODEL MinCVaR  'The CVaR model' /BudgetCon, ReturnCon, LossDef, VaRDevCon, ObjDefCVaR/;
MODEL MinCVaRLambda  'The Lambda CVaR model' /BudgetCon, ReturnCon, LossDef, VaRDevCon, ObjDefCVaR, ObjDefLambda/;
MODEL MinVaR  'The VaR model' /BudgetCon, ReturnCon, LossDef, VaRDevCon, SetZ, ObjDefVaR/;

scalar PortRet, WorstCase;
$ontext
MU_TARGET = MIN_MU;
SOLVE MinCVaR MINIMIZING CVaR USING LP;
PortRet = SUM(i, EP(i) * x.l(i));
VaR.l = -1*VaR.l;
CVaR.l = -1*CVaR.l;
WorstCase = -1*smax(l, Losses.l(l));
display VaR.l, CVaR.l, WorstCase, PortRet, x.l;

SOLVE MinVaR MINIMIZING VaR USING MIP;
PortRet = SUM(i, EP(i) * x.l(i));
VaR.l = -1*VaR.l;
WorstCase = -1*smax(l, Losses.l(l));
display VaR.l, WorstCase, PortRet, x.l;
*----------------------------------*


MU_TARGET = MAX_MU;
SOLVE MinCVaR MINIMIZING CVaR USING LP;
PortRet = SUM(i, EP(i) * x.l(i));
VaR.l = -1*VaR.l;
CVaR.l = -1*CVaR.l;
WorstCase = -1*smax(l, Losses.l(l));
display VaR.l, CVaR.l, WorstCase, PortRet, x.l;

SOLVE MinVaR MINIMIZING VaR USING MIP;
PortRet = SUM(i, EP(i) * x.l(i));
VaR.l = -1*VaR.l;
WorstCase = -1*smax(l, Losses.l(l));
display VaR.l, WorstCase, PortRet, x.l;
$offtext
*----------------------------------*



MU_TARGET = (MAX_MU+MIN_MU)/2;


display MU_TARGET;
MU_TARGET = MU_TARGET + 0.003;
display MU_TARGET;
lambda = 0.2;



SOLVE MinCVaRLambda MINIMIZING CVARLambda USING LP;
PortRet = SUM(i, EP(i) * x.l(i));
VaR.l = -1*VaR.l;
CVaR.l = -1*CVaR.l;
WorstCase = -1*smax(l, Losses.l(l));
display VarDev.l, VaR.l, CVaR.l, WorstCase, PortRet, x.l;




*SOLVE MinCVaR MINIMIZING CVaR USING LP;
*PortRet = SUM(i, EP(i) * x.l(i));
*VaR.l = -1*VaR.l;
*CVaR.l = -1*CVaR.l;
*WorstCase = -1*smax(l, Losses.l(l));
*display VarDev.l, VaR.l, CVaR.l, WorstCase, PortRet, x.l;
*
*
*
*
*
*SOLVE MinVaR MINIMIZING VaR USING MIP;
*PortRet = SUM(i, EP(i) * x.l(i));
*VaR.l = -1*VaR.l;
*WorstCase = -1*smax(l, Losses.l(l));
*display VarDev.l, VaR.l, WorstCase, PortRet, x.l;





SET FrontierPoints / PP_0 * PP_10 /

ALIAS(FrontierPoints,j);

PARAMETERS
         RiskWeight(j)           Investor's risk attitude parameter
         PortfolioCVAR(j)      Optimal level of portfolio variance
         PortfolioReturn(j)      Portfolio return
         OptimalAllocation(j,i)  Optimal asset allocation
         SolverStatus(j,*)       Status of the solver
         SummaryReport(*,*)      Summary report;

* The risk weight, lambda,  has to range in the interval [0,1]

RiskWeight(j) = (ORD(j)-1)/(CARD(j)-1);

StepSize = ((MAX_MU - MIN_MU) / (CARD(j)-1));
MU_TARGET = MIN_MU;
lambda = 0.2;
LOOP(j$(ord(j)>1 and ord(j)<card(j)),
   lambda = RiskWeight(j);
   MU_TARGET = MU_TARGET + StepSize;
   
    SOLVE MinCVaRLambda MINIMIZING CVARLambda USING LP;
    PortRet = SUM(i, EP(i) * x.l(i));
    VaR.l = -1*VaR.l;
    CVaR.l = -1*CVaR.l;
    WorstCase = -1*smax(l, Losses.l(l));
    display VarDev.l, VaR.l, CVaR.l, WorstCase, PortRet, x.l;

*   MeanVar.SOLPRINT = 2;
   PortfolioCVAR(j)= CVaR.l;
   PortfolioReturn(j)  = PortRet;
*   OptimalAllocation(j,i)     = x.l(i);
);

* Store results by rows

*SummaryReport(i,j) = OptimalAllocation(j,i);
SummaryReport('PortfolioCVAR',j) = PortfolioCVAR(j);
SummaryReport('PortfolioReturn',j) = PortfolioReturn(j);
SummaryReport('Lambda',j)   = RiskWeight(j);

DISPLAY SummaryReport;


EXECUTE_UNLOAD 'Summary.gdx', SummaryReport;
EXECUTE 'gdxxrw.exe Summary.gdx O=MeanCVARFrontier.xls par=SummaryReport rng=sheet1!a1' ;

