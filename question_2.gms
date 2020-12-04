$eolcom //
option optcr=0, reslim=120;

SET
         Scenarios
         Asset
;
ALIAS(Asset,i);
ALIAS(Scenarios, s);

PARAMETER
         ScenarioReturn(i, s)
;

$GDXIN FixedScenarios
$LOAD Asset, Scenarios, ScenarioReturn
$GDXIN


SCALARS
        Budget        'Nominal investment budget'
        alpha         'Confidence level'
        Lambda        'Risk aversion parameter'
        CVaRLim       'CVaR limit'
        ExpRetLim     'Expected Return Limit'
;

Budget = 100;
alpha  = 0.99;
CVaRLim = Budget*10;
ExpRetLim = -100;
lambda = 0.999;


PARAMETERS
        pr(s)       'Scenario probability'
        P(i,s)      'Final values'
        EP(i)       'Expected final values'
;

pr(s) = 1.0 / CARD(s);


P(i,s) = 1 + ScenarioReturn( i, s );


EP(i) = SUM(s, pr(s) * P(i,s));


POSITIVE VARIABLES
         x(i)            'Holdings of assets in monetary units (not proportions)'
         VaRDev(s)       'Measures of the deviation from the VaR'
;

VARIABLES
         losses(s)       'The scenario loss function'
         VaR             'The alpha Value-at-Risk'
         CVaR            'The alpha Conditional Value-at-Risk'
         ExpectedReturn  'Expected return of the portfolio'
         obj             'objective function value'
;

EQUATIONS
         BudgetCon       'Equation defining the budget constraint'
         ReturnCon       'Equation defining the portfolio expected return'
         LossDefCon(s)   'Equation defining the losses'
         VaRDevCon(s)    'Equation defining the VaRDev variable'
         CVaRDefCon      'Equation defining the CVaR'
         ObjectivFunc    'lambda formulation of the MeanCVaR model'
         CVaRLimCon      'Constraint limiting the CVaR'
         ReturnLimCon    'Constraint on a minimum expected return'
;



*--Objective------

*--s.t.-----------
BudgetCon ..             sum(i, x(i)) =E= Budget;

ReturnCon ..             ExpectedReturn =E= sum(i, EP(i)*x(i));

LossDefCon(s) ..         Losses(s) =E= -1*sum(i, P(i, s)*x(i) );

VaRDevCon(s) ..          VaRDev(s) =G= Losses(s) - VaR;

CVaRDefCon ..            CVaR =E= VaR + (sum(s, pr(s)*VarDev(s) ) )/(1 - alpha);

ObjectivFunc ..          Obj =E= (1-lambda)*ExpectedReturn - lambda*CVaR;

CVaRLimCon ..            CVaR =L= CVaRLim;

ReturnLimCon ..          ExpectedReturn =G= ExpRetLim;



*--Models-----------





MODEL CVaRModel 'The Conditional Value-at-Risk Model' /BudgetCon, ReturnCon, LossDefCon, VaRDevCon,CVaRDefCon, ObjectivFunc, CVaRLimCon, ReturnLimCon/;


*------------------Eff Front with equidistant CVaR--------------------*
set pp /pp0*pp10/;

PARAMETERS
         RunningCVaR(pp)          'Optimal level of portfolio CVaR'
         RunningReturn(pp)        'Portfolio return'
         RunningAllocation(pp,i)  'Optimal asset allocation'
         MaxCVar
         MinCVar
         Losses_all_scenarios(pp, s)
;

*first we find the biggest possible average, hence maximum CVaR:
lambda = 0.001;
SOLVE CVaRModel Maximizing OBJ Using LP;
MaxCVar = CVaR.l;
RunningCVaR(pp)$(ord(pp)=1) = CVaR.l;
RunningReturn(pp)$(ord(pp)=1)  = ExpectedReturn.l;
RunningAllocation(pp,i)$(ord(pp)=1)     = x.l(i);

*Then we find the lowest possible variance:
lambda = 0.999;
SOLVE CVaRModel Maximizing OBJ Using LP;
MinCVar = CVaR.l;
RunningCVaR(pp)$(ord(pp)=card(pp)) = CVaR.l;
RunningReturn(pp)$(ord(pp)=card(pp))  = ExpectedReturn.l;
RunningAllocation(pp,i)$(ord(pp)=card(pp))     = x.l(i);

display MaxCVar, MinCVar;

*Then we find the equidistant variances in between
lambda = 0.001;
CVarLim = MaxCVar;
loop(pp$(ord(pp)>1 and ord(pp)<card(pp)),
    CVarLim = CVarLim - (MaxCVar-MinCVar)/(card(pp)-1);
    SOLVE CVaRModel Maximizing OBJ Using LP;
*display VarLim, PortVariance.l;
    Losses_all_scenarios(pp,s) = Losses.l(s);
*    display Losses.l;
    RunningCVaR(pp) = CVaR.l;
    RunningReturn(pp)  = ExpectedReturn.l;
    RunningAllocation(pp,i)     = x.l(i);
);
display Losses_all_scenarios;
EXECUTE_UNLOAD 'Ques2_losses_all_scenarios.gdx', Losses_all_scenarios;
EXECUTE 'gdxxrw.exe Ques2.gdx O=Ques2_losses_all_scenarios.xls par=Losses_all_scenarios rng=sheet1!a1' ;

$exit
display RunningCVaR, RunningReturn, RunningAllocation;


parameter SummaryReport2(*,*);
* Store results by rows
SummaryReport2(i,pp) = RunningAllocation(pp,i);
SummaryReport2('CVaR',pp) = RunningCVaR(pp);
SummaryReport2('Return',pp) = RunningReturn(pp);


DISPLAY SummaryReport2;
* Write SummaryReport into an Excel file

EXECUTE_UNLOAD 'Ques2.gdx', SummaryReport2;
EXECUTE 'gdxxrw.exe Ques2.gdx O=MeanCVaRFrontier.xls par=SummaryReport2 rng=sheet1!a1' ;










//We fix the SPY's weight to the whole budget to obtain the CVaR and Target mean of SPY so that
//we can constrain our model accoring to SPY benckmark
X.fx(i) = 0;
X.fx('SPY') = Budget;
display X.l;

*------------CVaR----------------------
//
lambda = 0.999;
SOLVE CVaRModel Maximizing OBJ Using LP;


DISPLAY X.l, ExpectedReturn.l, VaR.L, CVaR.L;


Parameters
    MuTarget_SPY,
    CVaRTarget_SPY,
    BestCase,
    WorstCase
;

//Now we use the expected return and the CVaR of the equal weight strategy as target benchmarks in the CVaR model
MuTarget_SPY = ExpectedReturn.l;
CVaRTarget_SPY = CVaR.L;
display MuTarget_SPY, CVaRTarget_SPY;





//free the X variable 
X.lo(i) = 0;
X.up(i) = Budget*10;

lambda = 0.001;
CVaRLim =CVaRTarget_SPY;
ExpRetLim  = -100;
SOLVE CVaRModel Maximizing OBJ Using LP;
// these are the values that we are going 
DISPLAY X.l, ExpectedReturn.l, VaR.L, CVaR.L;



//free the X variable 
X.lo(i) = 0;
X.up(i) = Budget*10;
lambda = 0.999;
CVaRLim = Budget*10;
ExpRetLim = MuTarget_SPY;

SOLVE CVaRModel Maximizing OBJ Using LP;
// these are the values that we are going 
DISPLAY X.l, ExpectedReturn.l, VaR.L, CVaR.L;






