/*******************************************************************************

     PROGRAM:                              MLE Computer for NDMMF
                                           with Left Censoring

     Version:                              0.97

     SYSTEM:                               EMMMA

     AUTHOR:                               D. Donato

     DATE FIRST USED:                      February 2, 2006

     MODIFICATION DATES:                   October 5, 2010
                                           October 6, 2010
                                           October 7, 2010
                                           October 8, 2010
                                           October 9, 2010
                                           October 10, 2010
                                           October 11, 2010
                                           January 14, 2011
                                           January 17, 2011
                                           March 3, 2011
                                           April 13, 2011
                                           April 20, 2011
                                           May 16, 2011 - June 14, 2011
                                           June 24, 2011

     DESCRIPTION:

          This program computes the maximum likelihood estimators (MLE) of
          the numerous parameters of the National Model of Mercury in
          Fish Tissue (NDMMF) based on input data constituting
          a collection of uncensored and left-censored data. This program
          employs the Newton-Raphson method to determine best-fit parameters
          in the maximum-likelihood sense by iteration. Please refer to the
          external documentation for a full description of the equations
          and algorithm used in this program.

          This program is an adaptation of an MLE computer for a linear
          model. As in the case of the linear model, one of the parameters
          estimated by this program is the standard deviation (and therefore
          the variance) of the error term (the residuals) for the NDMMF.

          Although it is in general desirable to promote code reusability
          because programmer time is scarce and very expensive compared to machine
          time, the need for processing speed in this program outweighs
          the general (and arguably speculative) need for code reusability.
          Accordingly, global variables have been used freely in this program
          in order to reduce function-call overhead by avoiding the pushing and
          popping of non-static variables onto and off of the stack. In the next
          version, however, static variables will be used when possible in order
          to improve the readability of the code. This program's use of global
          variables is certainly unusual and admittedly contrary to what is normally 
          considered good programming style.
          

     MODIFICATION DESCRIPTIONS:

     January 17, 2011 -- The Newton-Raphson procedure for finding sigma squared
         was corrected.

     March 3, 2011 -- The version of this program successfully tested on
         January 17, 2011, assumed that censored observations were interval-
         censored with the detection limit as the upper bound on concentration
         and zero (0) as the lower bound -- this despite the identification
         of this program at the very top of this file as the "MLE Computer
         for NDMMF with Left Censoring". The mathematical expressions used
         in computing maximum-likelihood estimates for NDMMF parameters
         with interval censoring include a number of integrals of the
         normal PDF from a lower limit to an upper limit; for purposes
         of computation these integrals can be evaluated by finding
         the normal CDF at the upper limit and subtracting the value
         of the normal CDF at the lower limit. Although the computations
         are somewhat more complicated than this note suggests, it
         is essentially true that the assumption of left censoring
         means that the CDF is evaluated at a lower limit of negative
         infinity and is therefore equal to zero (0). Thus adapting prior
         code in this program for left censoring means dropping out
         some additive terms because they are now known in advance to
         be zero. Ideally this program should handle left, right, or
         interval censoring under the control of input directives or
         data, but this version has been specialized to left censoring
         in order to allow project work to move forward immediately.

    April 13, 2011 -- This version was modified to allow the use of an 
         integral weight for each observation. The weight must be greater
         than or equal to 1. If no observation has a weight other than 1,
         then there is, in effect, no weighting.

    April 20, 2011 -- The Newton-Raphson computations have been modified
         to perform multiple iterations for each parameter within a pass
         through all parameters to adjust the partial derivatives of
         the log-likelihood toward zero. A test has also been added
         within the Newton-Raphson code for each numerical estimation
         of a derivative to insure that the parameter is not changed
         if the derivative is very close to zero.

    May 16, 2011 - June 14, 2011 -- Various changes have been made in order
        to produce more precise computations of the standard normal CDF
        at extreme values, expecially extreme negative values, of x.

    June 24, 2011 -- The input of Hgdata.srt was modified to allow the
        input and storage of the integral record ID code.

*******************************************************************************/

/*******************************************************************************
                             PRE-PROCESSOR DIRECTIVES
*******************************************************************************/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define  SIZEOBSARRAY    120000 
#define  SIZEEVENTARRAY   20000
#define  SIZESPCARRAY       800
#define  NUMPARAMS        20000
#define  MAXITERS          2500


/*******************************************************************************
                                  GLOBAL VARIABLES
*******************************************************************************/


  /*** Mathematical Constants ***/

  double     zero={0.0}, one={1.0}, two={2.0}, four={4.0}, pi={3.141592654},
             pisqrt, twopisqrt, sqrtpidiv2, twosqrt, SmallTestValue={1.0e-12};

  /*** Data and Distribution Arrays and Related Variables ***/

  double     EPest[SIZEEVENTARRAY], prevEPest[SIZEEVENTARRAY],
             SPest[SIZESPCARRAY], prevSPest[SIZESPCARRAY];
  double     EPBest[SIZEEVENTARRAY],SPBest[SIZESPCARRAY];
  double     EPsave[SIZEEVENTARRAY], SPsave[SIZESPCARRAY];
  double     Resultin, Inchesin, Result[SIZEOBSARRAY], Lengthf[SIZEOBSARRAY];
  double     Paramestin, prevSigma, sigmain, numuncsr, wnumuncsr;
  int        SPC[SIZEOBSARRAY], Event[SIZEOBSARRAY], Dindex[SIZEOBSARRAY];
  int        RecordIDin, RecordID[SIZEOBSARRAY];
  double     WT[SIZEOBSARRAY], WTin;
  int        SPCin, Eventin, DLin, Dindexin;
  int        obscount, eventcount, spccount;
  unsigned char  DL[SIZEOBSARRAY], EPindic[SIZEEVENTARRAY], SPindic[SIZESPCARRAY];

  int        cOBS, cSPC, cEvent;  /* Current values of these indexes. */
  int        cOBSnext, cSPCnext, cEventnext; /* Next values of these indexes. */
  int        rEvent, rSPC;        /* Loop variables for events and SPC's,
                                     respectively. */

  /*** Variables for Function Return Values ***/

  double      ndpdf, ndcdf, dwrta, dwrtb, dwrts2, dFwrts2, dlikelihood,
              dloglikelihood, altA, altAcount;

  /*** Variables for computing the Standard Normal CDF. ***/

  double      xsn, sndcdf, sncumdev[800001];
  int         extremecount;

  /*** Variables for approximating dwrta and dwrtb. ***/

  int        sigmagridlevel;
  double     xhigh, xlow, xrange, x0, xnc, d4, d5;

  /*** Variables for computing dFwrts2.             ***/

  double     AdFwrts2, ds2cum[800001], ds2cummax;

  /*** Variables for Computing the Sample Variance of the Error Term ***/

  double     xsum, x2sum, xcount;

  /*** Utility Variables and Throw-Away Intermediate Values ***/

  double     temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, incrementA;
  double     temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18;
  double     temp21, temp22, temp23, temp24, temp25, temp26, temp27, temp28;
  int        i, j, k, m, ia, ib, ic, il, inr, inrnext, jnr, retval;

  /*** General Computational Variables ***/

  double      x, mu, a, b, g, sigma, sigma2, sigmaBest;
  int         intnum;

  /*** Globalized Function Variables ***/

  /** ndevcdf **/
  double xt;

  /** derivwrts2, likelihood, and loglikelihood **/
  double sum0, sum1, sum2, prod1;

  /** numerical derivatives of derivatives **/
  double numdrva, numdrvb, numdrvs2;

  /** Newton-Raphson **/
  double nrtp1, nrtp2, nrtp3, nrtp4, nrtp5, nrtp6;
  double dwrtasum, dwrtbsum, dwrtacount, dwrtbcount, sigmasave;
  double ddiv, Ederivincrement[SIZEEVENTARRAY],
         Sderivincrement[SIZESPCARRAY], sigmaderivincrement;
  double maxlikelihood, maxloglikelihood, llsave[MAXITERS + 1];
  int maxliteration, maxlogiteration, tndx;
  int numiters, indicator, EPtype, SPtype;
  unsigned char algsign;

  /*** File-Related Variables ***/
  FILE  *OUT, *EP, *SP, *IN1, *IN2, *IN3, *IN4, *IN5, *ofp;

/*******************************************************************************
                                 FUNCTION PROTOTYPES
*******************************************************************************/


  /*** General Processing and I/O Functions ***/

  void      Initialization(void);
  void      DataInput(void);
  void      IdentifyExcludedEPs(void);
  void      IdentifyAllCensoredEPs(void);
  void      IdentifyExcludedSPs(void);
  void      ReportResults(void);
  void      WriteBestParams(void);


  /*** Mathematical and Computational Functions ***/

  void    ndevpdf(void);         /* Computes double from x, mu, sigma.        */
  void    ndevcdf(void);         /* Computes double from x, mu, sigma.        */
  double  sndevcdfapprox(double);/* Computes double from xt.                  */
  double sndevcdfapproxB(double);/* Computes double from xt.                  */
  double sndevcdfapproxC(double);/* Computes double from xt.                  */
  void    sndevcdfinit(void);    /* Initializes Normal CDF computation.       */
  double  approxarea(double, double);
  double  approxareaB(double, double);
  void    sndevcdf(void);        /* Computes double from fabs((x-mu)/sigma).  */
  void    derivwrta(void);       /* Computes double from x.                   */
  void    approxdwrtab(void);    /* Computes double from xhigh, xlow, etc.    */
  void    computewrtab(void);    /* Computes double from xhigh and xlow.      */
  void    derivwrtb(void);       /* Computes double from b.                   */
  void    derivwrts2(void);      /* Computes double from sigma squared.       */
  void    derivFwrts2(void);     /* Computes double from x.                   */
  void    AderivFwrts2(void);    /* Computes double from x.                   */
  void    ds2cuminit(void);      /* Computes numerical integral for dFwrts2.  */
  void    alternativeA(void);    /* Computes double from s, sigma;            */
  void    numderiva(void);       /* Computes double from cOBS.                */
  void    numderivb(void);       /* Computes double from cSPC.                */
  void    numderivs2(void);      /* Computes double from sigma2               */
  void    likelihood(void);      /* Computes double from a, b, and sigma.     */
  void    loglikelihood(void);   /* Computes double from a, b, and sigma.     */
  void    SigmaGrid(double);     /* Computes best sigma using a grid.         */
  void    SigmaSign(double);     /* Computes best sigma by sign change.       */
  void    NewtonRaphson(void);   /* Computes parameter estimates by iterative */
                                 /*    successive approximation.              */

  /*** Test Function(s)                          ***/

  void    CompDenomTerm(void);  /* Compute the problematic denominator term. */

/*******************************************************************************
                               EXECUTIVE (MAIN) ROUTINE
*******************************************************************************/

int  main(int argc, char * argv[])
{
printf("Starting data input.\n");
  DataInput();
printf("Starting initialization. \n");
  Initialization();
/*** Temporary Call        ****
  CompDenomTerm();
**** End of Temporary Call ***/
printf("Starting NewtonRaphson.\n");
  NewtonRaphson();
  ReportResults();
  WriteBestParams();

  return(0);
}

/*******************************************************************************
                                INITIALIZATION ROUTINE
*******************************************************************************/

  void  Initialization()
{
  IdentifyAllCensoredEPs();
  IdentifyExcludedEPs();
  IdentifyExcludedSPs();

  pi = M_PI;
  twosqrt=sqrt(2.0);
  pisqrt=sqrt(pi);
  sqrtpidiv2=sqrt(pi/2.0);
  twopisqrt=sqrt(two*pi);
  mu=0.0;
  sndevcdfinit();
  ds2cuminit();

  sigma=sigmain; sigma2=sigma*sigma; sigmaBest=sigma;

  maxlikelihood=0.0; maxliteration=0;
  maxloglikelihood=-9.9E10; maxlogiteration=0;
  sigmaderivincrement=0.001;

  printf("eventcount and spccount are: %5d  %5d\n", eventcount, spccount);
  printf("numuncsr is: %5lf\n", numuncsr);

  return;
}

/*******************************************************************************
                                  DATA INPUT ROUTINE
*******************************************************************************/

  void DataInput()
{
  double tot, tot2;

  /** Zero out the input arrays. **/

  for (i=0; i<SIZEOBSARRAY; i++)
       {SPC[i]=0; Event[i]=0; Result[i]=0.0; Lengthf[i]=0.0; DL[i]=0; Dindex[i]=0;}
  for (i=0; i<SIZEEVENTARRAY; i++)
       {EPest[i]=0.0; prevEPest[i]=0.0; EPsave[i]=0.0; EPindic[i]=0;}
  for (i=0; i<SIZESPCARRAY; i++)
       {SPest[i]=0.0; prevSPest[i]=0.0; SPsave[i]=0.0; SPindic[i]=0;}

  /** Input observations sorted by event and the inversion pointer used  **/
  /** to provide access to observations by SPC.                           **/

  obscount=0; numuncsr=0.0; wnumuncsr=0.0;
  IN1=fopen("Hgdata.srt","r");
  IN2=fopen("Hgdata.ndx","r");
  tot=0.0; tot2=0.0;
  while(!feof(IN1))
    {
     retval=fscanf(IN1,"%d %d %lf %lf %d %lf %d\n",
                   &SPCin, &Eventin, &Inchesin, &Resultin, &DLin, &WTin, &RecordIDin);
     if (retval!=7) {break;}
     if (DLin==0) {numuncsr=numuncsr+1.0; wnumuncsr+=(double)WTin;}
     obscount++;
     retval=fscanf(IN2,"%d\n", &Dindexin);
     SPC[obscount]=SPCin; Event[obscount]=Eventin; DL[obscount]=DLin;
     WT[obscount] = WTin; RecordID[obscount] = RecordIDin;
     Result[obscount]= log((Resultin*1000.0) + 1.0);
     tot+=Result[obscount]; tot2+=Result[obscount]*Result[obscount];
     Lengthf[obscount]=log(Inchesin + 1.0);
     Dindex[obscount]=Dindexin;
    }

  fclose(IN1); fclose(IN2);

printf("tot, tot2, N, and avg: %12.4lf | %12.4lf | %5d | %12.8lf\n",
        tot, tot2, obscount, tot/(double)obscount);
printf("Sample SD: %12.8lf\n",sqrt( ( tot2-(tot*tot/(double)obscount) ) / (double) (obscount-1)));
printf("Weighted number of uncensored observations is: %12.2lf\n", wnumuncsr);

  /** Input initial estimates for event parameters. **/

  eventcount=0;
  IN3=fopen("Hgevents.srt","r");

  while(!feof(IN3))
    {
     retval=fscanf(IN3,"%d %lf\n", &Eventin, &Paramestin);
     if (retval!=2) {break;}
     EPest[Eventin]=Paramestin;
     EPBest[Eventin]=Paramestin;
     EPsave[Eventin]=Paramestin; EPindic[Eventin]=1;
     prevEPest[Eventin]=Paramestin-0.017528;
     Ederivincrement[Eventin]=0.001;
     eventcount++;
    }

  fclose(IN3);

  /** Input initial estimates for species/cut (SPC) parameters. **/

  spccount=0;
  IN4=fopen("Hgspc.srt","r");

  while(!feof(IN4))
    {
     retval=fscanf(IN4,"%d %lf\n", &SPCin, &Paramestin);
     if (retval!=2) {break;}
     SPest[SPCin]=Paramestin;
     SPBest[SPCin]=Paramestin;
     SPsave[SPCin]=Paramestin; SPindic[SPCin]=1;
     prevSPest[SPCin]=Paramestin-0.00001;
     Sderivincrement[SPCin]=0.001;
     spccount++;
    }

  fclose(IN4);

  /** Input initial estimate for sigma. **/

  IN5=fopen("Hgsigma.dat","r");
  fscanf(IN5, "%lf\n", &sigmain);
  fclose(IN5);

  return;
}

/*******************************************************************************
             ROUTINE FOR TEST COMPUTATION OF THE PROBLEMATIC DENOMINATOR
*******************************************************************************/

  void  CompDenomTerm()
{
 int i;
 double x1, increment;

 x1 = 10.0;
 increment = 0.1;

 while (x1 >= -15.0)
  {
   x = x1+four;
     ndevcdf(); temp1=ndcdf;
     ndevpdf(); temp2=ndpdf;
   x = x1;
     ndevcdf(); temp3=ndcdf;
     ndevpdf(); temp4=ndpdf;
     temp5=temp4-temp2; temp6=temp1-temp3;
     temp7=log(fabs(temp5))-log(fabs(temp6));
     temp7=exp(temp7); if (temp5<zero) {temp7=-temp7;}
   printf("<><><>%16.13lf %16.13lf %16.13lf %16.13lf %16.13lf \n", x1, x1+one, temp5, temp6, temp7);
   x1-=increment;
  }

 return;
}
/*******************************************************************************
          IdentifyExcludedEPs: Find the events for which data exist in the
               data input file but not in the event-parameter input file
               and mark them for exclusion by setting EPindic[cEvent] to 5.
*******************************************************************************/

  void IdentifyExcludedEPs()
{
  for (cOBS=1; cOBS<=obscount; cOBS++)
      {
       cEvent=Event[cOBS];
       if (EPindic[cEvent]==0) {EPindic[cEvent]=5;}
/*if (DL[cOBS]!=0) {EPindic[cEvent]=5;} */
      }
  return;
}

/*******************************************************************************
          IdentifyAllCensoredEPs: Find the events for which all observations
               are left censored and mark them for simple algebraic
               computation of their respective Event Parametersi by setting
               EPindic[cEvent] to 8.
*******************************************************************************/

  void IdentifyAllCensoredEPs()
{
  int prevcEvent, allcount, csrcount;

  for (cOBS=1; cOBS<=obscount; cOBS++)
      {
       cEvent=Event[cOBS];
       if (cOBS==1) {prevcEvent=cEvent; allcount=0; csrcount=0;}
       if (cEvent!=prevcEvent)
            {
             if (EPindic[prevcEvent]==1 && allcount==csrcount)
                  {EPindic[prevcEvent]=8;}
             prevcEvent=cEvent; allcount=0; csrcount=0;
            }
       allcount++; if (DL[cOBS]==1) {csrcount++;}
      }
  if (EPindic[prevcEvent]==1 && allcount==csrcount)
      {EPindic[prevcEvent]=8;}

  return;
}

/*******************************************************************************
          IdentifyExcludedSPs: Find the SPCs for which data exist in the
               data input file but not in the SPC-parameter input file
               and mark them for exclusion by setting SPindic[cSPC] to 5.
*******************************************************************************/

  void IdentifyExcludedSPs()
{
  for (cOBS=1; cOBS<=obscount; cOBS++)
      {
       cSPC=SPC[cOBS];
       if (SPindic[cSPC]==0) {SPindic[cSPC]=5; printf("***Excluded SPC: %4d\n", cSPC);}
      }
  return;
}

/*******************************************************************************
                          RESULT REPORTING ROUTINE
*******************************************************************************/

  void ReportResults()
{
  printf("\n\n**************************************************************\n\n");
  printf("The event parameter values determined by MLE are:\n\n");
  for (i=0; i<SIZEEVENTARRAY; i++)
    {
     if (EPindic[i]==1 || EPindic[i]==8)
     printf("[%5d] ==> %10.6lf %10.6lf | %12.9lf\n",i, EPsave[i], EPBest[i], EPsave[i]-EPBest[i]);
     if (EPindic[i]==5)
     printf("[%5d] *** Excluded from computations!\n", i);
    }
  printf("\n\n**********************************************\n\n");

  printf("The species/cut (SPC) parameter values determined by MLE are:\n\n");
  for (i=0; i<SIZESPCARRAY; i++)
    {
     if (SPindic[i]==1 || SPindic[i]==8)
     printf("[%5d] ==> %10.6lf %10.6lf |  %12.9lf\n",i, SPsave[i], SPBest[i], SPsave[i]-SPBest[i]);
     if (SPindic[i]==5)
     printf("[%5d] *** Excluded from computations!\n", i);
    }
  printf("\n\n**********************************************\n\n");

  printf("The  estimate of the standard-deviation parameter determined by MLE is:\n\n");
     printf("[-----] ==>   %14.10lf\n", sigmaBest);
  printf("\n\n**********************************************\n\n");

  printf("The computation used %3d iterations.\n\n", numiters);

  printf ("The maximum loglikelihood computed was %14.12lf on iteration number %4d.\n\n",
          maxloglikelihood, maxlogiteration);
    /*MNF added file with final results*/
    ofp = fopen("summaryRESULTS.dat","w");
    fprintf (ofp, "%20s %20s %20s %20s\n" , "max sig","max loglike","total iters","best iteration");
    fprintf (ofp, "%20f %20f %20d %20d\n", sigmaBest, maxloglikelihood,numiters,maxlogiteration);
    fclose(ofp);
  return;
}

/*******************************************************************************
                          WRITE THE BEST PARAMETERS TO FILES
*******************************************************************************/

void  WriteBestParams()
{
 EP = fopen("BestEPs", "w");
 for (i=0; i<SIZEEVENTARRAY; i++)
    {
     if (EPindic[i] == 1 || EPindic[i] == 8)
     fprintf(EP, " %5d  %12.9lf\n", i, EPBest[i]);
    }
 fclose(EP);

 SP = fopen("BestSPs", "w");
 for (i=0; i<SIZESPCARRAY; i++)
    {
     if (SPindic[i]==1 || SPindic[i]==8)
     fprintf(SP, " %5d  %12.9lf\n",i, SPBest[i]);
    }
 fclose(SP);

 return;
}

/*******************************************************************************
          ndevpdf :  Compute a Value of the Normal PDF.
*******************************************************************************/

  void ndevpdf()                    /* Computes double from x, mu, and sigma. */
{
  ndpdf=exp(-(pow((x-mu),two)/(two*sigma2))) / (sigma*twopisqrt);
  return;
}

/*******************************************************************************
          ndevcdf : Compute a value of the Normal CDF.
*******************************************************************************/

  void ndevcdf()                    /* Computes double from x, mu, and sigma. */
{
  static long double xtl, xl, mul, sigmal, ndcdfl, twol, point5l;
  twol = 2.0; point5l = 0.5;
  xl = x;
  mul = mu;
  sigmal = sigma;

  m=0;
/*xt=(x-mu)/sigma; */
  xtl = (xl - mul)/sigmal;

/*if (xt<zero) {xt=-xt; m=1;} */
  if (xt<zero) {xtl=-xtl; m=1;}

  if (xt <= 16.00000 || (m==0) )
    {
     ndcdfl=erfl(xtl/sqrtl(twol))/twol;
     if (m==0) {ndcdfl=point5l+ndcdfl;}
     else {ndcdfl=point5l-ndcdfl;}
     ndcdf=ndcdfl;
    }
  else
    {
     xt = -xt;
     sndevcdf();
     ndcdf = xsn;
    }

  return;
}

/*******************************************************************************
          sndevcdfapprox : Compute an approximation of the value of the standard
                           normal CDF at xt.
*******************************************************************************/

  double sndevcdfapprox(double xt)
{
  double approx, cumvalue, xsave, xt2;
  int    altsign, ii;

  xsave = x;

  xt2 = xt * xt;
  approx = 1.0;
  cumvalue = 0.0;
  altsign = -1;

  for (ii=0; ii<30; ii++)
     {
      cumvalue -= log(xt2);   
      cumvalue += log((double)((ii*2)+1));
      if (altsign < 0)   {approx -= exp(cumvalue);}
      else               {approx += exp(cumvalue);}
      if (altsign == -1) {altsign = 1;} else {altsign = -1;}
     }

  ndevpdf();
  approx*= (ndpdf / xt);

  x = xsave;
  return(approx);
}

/*******************************************************************************
        sndevcdfapproxB : Compute a power-series  approximation of the value
                          of the standard normal CDF at xt.
*******************************************************************************/

  double sndevcdfapproxB(double xt)
{
  double approx, iid, prod, sum, xsave, xt2;
  int    altsign, ii;

  xt2 = xt * xt;
  prod  = xt;
  sum   = xt;
  altsign = -1;

  for (ii=1; ii<10000; ii++)
     {
      iid = ii;
      prod *= xt2 *(two*(iid-one) + one) / ( iid * two * (two*iid + one));
/*    if (fabs(prod) < 1.0E-90) {printf("%5d\n", ii); break;} */
      if (altsign < 0)   {sum -= prod;}
      else               {sum += prod;}
      if (altsign == -1) {altsign = 1;} else {altsign = -1;}
printf("ii prod sum: %5d %12.9lf %12.9lf\n", ii, prod, sum);
     }

  approx = (one/two) + (one / twopisqrt) * sum;
  return(approx);
}

/*******************************************************************************
         sndevcdfapproxC: Compute an approximation of the value of the standard
                          normal CDF at xt using Marsaglia's method.
*******************************************************************************/

  double sndevcdfapproxC(double xt)
{
  double approx, cumvalue, xsave, xt2;
  int    ii;

  xsave = x;

  xt  = xt + 6.0;
  xt2 = xt * xt;
  approx = xt;
  cumvalue = xt; 

  for (ii=1; ii<300000; ii++)
     {
      cumvalue *= (xt2 / ( (double) ((ii * 2) + 1)  ));   
      approx += cumvalue;
/*printf("xt approx cumvalue: %12.9lf %18.14E  %18.14E\n", xt, approx, cumvalue);*/
     }

  x = xt -6.0;
  ndevpdf();
  approx *= ndpdf;
  approx += (double) 0.5 + erf(6.0/sqrt(2.0))/2.0;

  x = xsave;
  return(approx);
}


/*******************************************************************************
          sndevcdfinit : Initialize Standard Normal CDF computation.
*******************************************************************************/

/** N.B. This function has been modified to compute the tail of the CDF for  **/
/**      extreme negative values of x between -8 and -16 standard units      **/
/**      below the mean.    June 5, 2011 -- D. Donato                        **/

  void sndevcdfinit()
{
  double v1, v2;

  sncumdev[0] = 0.5 - (erfl(8.0/sqrt(2.0))/2.0);
  x = 8.0; mu = 0.0;
  sigma = 1.0; sigma2 = 1.0;
  ndevpdf(); v1 = ndpdf;
  x = 8.0001;

  for (i=1; i<=80000; i++)
     {
      ndevpdf(); v2 = ndpdf;
/*    sncumdev[i] = sncumdev[i-1] - (((v1+v2)/two)*0.00001); */
      sncumdev[i] = sncumdev[i-1] - approxareaB(x-0.0001, x);
/*printf("x v1 v2 v1-v2: %10.5lf %18.14E %18.14E %18.14E\n", x, v1, v2, v1-v2);*/
/*printf("i sncumdev[i]: %6d %18.14E\n", i, sncumdev[i]);  */
      x = 8.0 + ((double)(i+1)) * 0.0001;
      v1 = v2;
     }

  return;
}

/*******************************************************************************
          approxarea: Compute the approximate area under the PDF curve
                      for an interval.
*******************************************************************************/

  double approxarea(double x1, double x2)
{
 double a, b, c, s, xsave, y1, y2;

 xsave = x;
 x = x2; ndevpdf(); y2 = ndpdf;
 x = x1; ndevpdf(); y1 = ndpdf;
 s = -x * y1;
 x = xsave;

 a = ((y1-y2) - s * (x1-x2)) / ((x1*x1 - x2*x2) - two*x1*(x1-x2));
 b = s - (two * a * x1);
 c = y1 - (a*x1*x1) - (b*x1);

 return( ((a*x2*x2*x2/3.0) + (b*x2*x2/2.0) + c*x2 )
        -((a*x1*x1*x1/3.0) + (b*x1*x1/2.0) + c*x1 ) );
}


/*******************************************************************************
         approxareaB: Compute the approximate area under the PDF curve
                      for an interval.
*******************************************************************************/

  double approxareaB(double x1, double x2)
{
  double     A1, A2, A3, Ae, a2, f1, f2, b1, b2, w, xmid, xsave;
  int        i, j, k;

  xsave = x;
  w  = x2 - x1;

     x = x1; ndevpdf(); f1 = ndpdf;
     x = x2; ndevpdf(); f2 = ndpdf;
     b1 = (f2 - f1)/w; A1 = ((f2 + f1)/two)*w;
     xmid = x1 + ((x2-x1)/two); x = xmid; ndevpdf(); 
     b2 = -(xmid)*ndpdf; a2 = ndpdf - (b2 * xmid);
     A2 = (((a2 + b2 * x1) + (a2 + b2 * x2))/two) * w;
     A3 = (A1 + A2)/two;
     Ae = (A1 -A2)/two;

 x = xsave;
/*printf("A3 Ae: %18.14E %18.14E\n", A3, A3);*/
 return(A3);

}

/*******************************************************************************
          sndevcdf : Compute a value of the Standard Normal CDF.
*******************************************************************************/

/** N.B. This function has been modified for use only for values of xt       **/
/**      between -8.0 and -16.0 standard units.  June 5, 2011 -- D. Donato   **/

  void sndevcdf(void)        /* Computes double from fabs((x-mu)/sigma).  */
{
  double cuminc, xpos, xposinc, smallnum;
  int    index1;
  smallnum = exp(-1.0E12);

  xpos = fabs(xt + 8.0);
 
  index1  = xpos * 10000.0;
  cuminc  = sncumdev[index1 + 1] - sncumdev[index1];
  xposinc = xpos - (((double)(index1))/10000.0);
  xsn     = sncumdev[index1] + ((xposinc) * cuminc);
  return;
}

/*******************************************************************************
          derivwrta :  Compute the derivative of the log-likelihood function
                       with respect to (wrt) one of the (alpha) event
                       parameters.
*******************************************************************************/

  void derivwrta()
{
  sigma2=sigma*sigma;
  dwrta=zero;
  rEvent = cEvent;
  cOBSnext=cOBS;

  while (rEvent == cEvent)
    {
     x=Result[cOBSnext]-EPest[cEvent]-SPest[SPC[cOBSnext]]*Lengthf[cOBSnext];
     xsum += x; x2sum += (x*x); xcount++;  /* Accumulate values for estimating sigma. */

     if (DL[cOBSnext]==0) /* If the observation is uncensored... */
        {
         dwrta+=(x/sigma2)*((double)WT[cOBSnext]);
        }
     else                 /* If the observation is censored... */
        {
         xhigh = x;
         computewrtab();
         dwrta += temp7*((double)WT[cOBSnext]);
/*
         xhigh = x;
         xlow  = -EPest[cEvent] - SPest[SPC[cOBSnext]]*Lengthf[cOBSnext];

         ** Determine whether to compute dwrta or approximate it **
         ** based on how extreme the high and low x values are.  **
         xrange=xhigh-xlow; x0 = -xrange/two;
         if (xlow > x0 + ((double)7.0 * sigma) || xhigh < x0 - ((double)7.0 * sigma))
           {
            ** Choose to approximate. **
            approxdwrtab();
            dwrta += temp7*((double)WT[cOBSnext]);
printf("Approx: x temp7 WT -- %12.9lf %12.9lf %3d\n", x, temp7, WT[cOBSnext]);
           }
         else
           {
            ** Compute directly. **
            computewrtab();
            dwrta += temp7*((double)WT[cOBSnext]);
**printf("Exact:  x temp7 WT -- %12.9lf %12.9lf %3d\n", x, temp7, WT[cOBSnext]);**
           }
*/
        }
     cOBSnext++;
     if (cOBSnext<=obscount) { rEvent=Event[cOBSnext];}
     else {rEvent=-1;}
    }

  return;
}

/********************************************************************************
          approxdwrta: Compute an approximation for the contribution of a
                       censored observation to the derivative of the
                       log-likelihood function with respect to (wrt) one of
                       the (alpha) event parameters. This approximation is
                       only called when one or both of the high or low x values
                       for the observation are extreme (that is, a large
                       multiple of the standard deviation above or below 
                       the mean for the error term).
*******************************************************************************/

  void approxdwrtab()
{
  double temp1, temp2, temp3, temp4, temp5, temp6;

  /* xhigh, xlow, x0, xrange, and sigma must be set before calling this */
  /* function.                                                          */

  if (xlow > x0 + ((double)7.0 * sigma))
    {
     /** Approximate a value in the high positive non-computable area. **/
     xnc = xlow;
     xhigh = ((double)5.0*sigma + x0); xlow = xhigh - xrange;
     computewrtab(); d5 = temp7;
     xhigh = ((double)4.0*sigma + x0); xlow = xhigh - xrange;
     computewrtab(); d4 = temp7;
     temp7 = d5 + ((d5-d4)/sigma)*(xnc-((double)5.0*sigma + x0));
    }
  else
    {
     /** Approximate a value in the low negative non-computable area. **/
     xnc = xhigh;
     xlow = (-(double)5.0*sigma - x0); xhigh = xlow + xrange;
     computewrtab(); d5 = temp7;
     xlow = (-(double)4.0*sigma - x0); xhigh = xlow + xrange;
     computewrtab(); d4 = temp7;
     temp7 = d5 - ((d4 - d5)/sigma)*(xnc + ((double)5.0*sigma + x0));
    }
 
  return;
}

/*******************************************************************************
          computedwrtab: Compute the contribution of a censored observation to
                         the derivative of the log-likelihood with respect to
                         (wrt) one of the alpha or beta parameters. In the
                         case of alpha parameters, the value computed by this
                         function is the contribution; in the case of the beta
                         parameters, the value computed by this function must
                         be multiplied by "g" to produce the actual
                         contribution. The contribution value is returned
                         in temp7.
*******************************************************************************/

  void computewrtab()
{
  double temp1, temp2, temp3, temp4, temp5, temp6;

  /* xhigh and xlow must be set before calling this function. */

  x = xhigh;
    ndevcdf(); temp1=ndcdf;
    ndevpdf(); temp2=ndpdf;
    if (temp1 == 0.0) {temp1 = ndpdf / fabs(x);}
/*** March 3, 2011 -- Changed from interval to left censoring.
  x = xlow;
    ndevcdf(); temp3=ndcdf;
    ndevpdf(); temp4=ndpdf;
***/
/*March 3, 2011 - Left censor*/ temp3=(double)0.0; temp4=(double)0.0;
    temp5=temp4-temp2; temp6=temp1-temp3;
    temp7=log(fabs(temp5))-log(fabs(temp6));
    temp7=exp(temp7);
      if (temp5<zero) {temp7=-temp7;}
      if (temp6<zero) {temp7=-temp7;}
    if(isinf(temp7) || isnan(temp7))
      {temp7=-fabs(x); printf("Problem computing dwrta or dwrtb.\n");}
  return;
}

/*******************************************************************************
          derivwrtb :  Compute the derivative of the log-likelihood function
                       with respect to (wrt) one of the (beta) SPC parameters.
*******************************************************************************/

  void derivwrtb()
{
  sigma2=sigma*sigma;
  dwrtb=zero;
  rSPC = cSPC;
  inrnext=inr;

  while (rSPC == cSPC)
    {
     tndx=Dindex[inrnext];
     g=Lengthf[tndx]; 
     x=Result[tndx]-EPest[Event[tndx]]-SPest[cSPC]*g;
     if (DL[tndx]==0) /* If the observation is uncensored... */
        {
         dwrtb+=((g*x)/sigma2)*((double)WT[tndx]);
        }
     else                 /* If the observation is censored... */
        {
         xhigh = x;
         computewrtab();
         dwrtb += g*temp7*((double)WT[tndx]);
/*
         xlow  = -EPest[Event[tndx]] - SPest[cSPC]*g;
         ** Determine whether to compute dwrta or approximate it **
         ** based on how extreme the high and low x values are.  **
         xrange=xhigh-xlow; x0 = -xrange/two;
         if (xlow > x0 + ((double)7.0 * sigma) || xhigh < x0 - ((double)7.0 * sigma))
           {
            ** Choose to approximate. **
            approxdwrtab();
            dwrtb += g*temp7*((double)WT[tndx]);
           }
         else
           {
            ** Compute directly. **
            computewrtab();
            dwrtb += g*temp7*((double)WT[tndx]);
           }
*/
        }
     inrnext++;
     if (inrnext<=obscount) { rSPC=SPC[Dindex[inrnext]];}
     else {rSPC=-1;}
    }

  return;
}

/*******************************************************************************
          derivwrts2 : Compute the derivative of the log-likelihood function
                       with respect to (wrt) sigma squared (the variance of
                       the error term). Also preserve the summation terms.
*******************************************************************************/

  void derivwrts2()
{
  double     temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8;

  altAcount = 0.0;
  extremecount = 0;
  sigma2=sigma*sigma;
  sum1=0.0; sum2=0.0;              /* sum1 and sum2 are used later            */
  dwrts2=-wnumuncsr/(two*sigma2);  /* in the likelihood() and loglikelihood() */
  cOBS=1;                          /* functions.                              */

  while (cOBS<=obscount)
    {
     cEvent=Event[cOBS]; cSPC=SPC[cOBS];
     if (EPindic[cEvent]!=1 || SPindic[cSPC]!=1) {cOBS++; continue;}
     x=Result[cOBS]-EPest[cEvent]-SPest[cSPC]*Lengthf[cOBS];
if (fabs(x/sigma)  > 6.0) extremecount++;
/*if (fabs(x/sigma) > 5.0) printf("RecordID SPC Event x/sigma DL: %6d %4d %5d %12.9lf %d\n", RecordID[cOBS], SPC[cOBS], Event[cOBS], x/sigma, DL[cOBS]);*/
/*if (fabs(x/sigma) > 8.0) printf("x = %12.9lf\n", x/sigma); */
     if (DL[cOBS]==0)
          {
           if (!isinf(x) && !isnan(x)) {sum1+=(x*x)*((double)WT[cOBS]);}
           else {printf("Error: x is invalid.\n");}
          }
     else
          {
           ndevcdf(); temp1=ndcdf; derivFwrts2(); temp2=dFwrts2;
/* if (fabs(x/sigma) > 6.0) printf("cOBS x x/sigma ndcdf(x) log(temp1) = %5d %12.9lf %12.9lf %18.16e %12.9lf\n", cOBS, x, x/sigma, temp1, log(temp1)); */
/*printf("log(1.0e-50): %12.9lf\n", log(1.0e-50)); */
           AderivFwrts2();
           if (fabs(AdFwrts2 - temp2) > 0.0005)
             {
              temp2 = AdFwrts2;
              printf("dFwrts2 AdFwrts2: %12.9lf %12.9lf \n", dFwrts2, AdFwrts2);
             }
/*** March 3, 2011 - Changed from interval to left censoring.
           x=-EPest[cEvent]-SPest[cSPC]*Lengthf[cOBS];
           ndevcdf(); temp3=ndcdf; derivFwrts2(); temp4=dFwrts2;
***/
           if ((temp1 <= 0.0) || (isinf(temp1)) || (isinf(temp1))) {ndevpdf(); temp1 = ndpdf /fabs(x);}
           temp3=(double)0.0; temp4=(double)0.0;
           temp5=temp2-temp4; temp6=temp1-temp3;
           temp7=log(fabs(temp5))-log(fabs(temp6));
           temp7=exp(temp7);
             if (temp5<zero) {temp7=-temp7;}
             if (temp6<zero) {temp7=-temp7;}
           if (!isinf(temp7) && !isnan(temp7))
             {
              dwrts2 += temp7*((double)WT[cOBS]);
             }
           else
             {
              alternativeA(); dwrts2 += altA * ((double)WT[cOBS]);
              altAcount++;
             }
           if (isnan(log(temp1))) temp1 = 1.0E-20; 
           sum2 += log(temp1) * ((double)WT[cOBS]);
          }
     cOBS++;
    }

    sum1/=two*sigma2*sigma2;
    dwrts2+=sum1;

  return;
}

/*******************************************************************************
          derivFwrts2 : Compute the derivative of the CDF with respect to
                        sigma squared at x.
*******************************************************************************/

  void derivFwrts2()
{
  double temp11, temp12, temp13, incrementA;
  int    ii;

  sigma2=sigma*sigma;
  temp11 = sigma;

  for (ii=0; ii<3; ii++)
    {
     incrementA= 0.001 + ((double)ii) * 0.005;
     sigma2 += incrementA; sigma=sqrt(sigma2); ndevcdf(); temp12=ndcdf;
     sigma2 -= two*incrementA; sigma=sqrt(sigma2); ndevcdf(); temp13=ndcdf;
     dFwrts2 = (temp12-temp13)/(two*incrementA);
     if (dFwrts2 != 0.0) {break;}
    }

  sigma = temp11; sigma2=sigma*sigma;

/*AderivFwrts2();
if (fabs(x) > 4.0) printf("x, sigma, dFwrts2, AdFwrts2: %12.9lf %12.9lf %14.10E %14.10E \n", x, sigma, dFwrts2, AdFwrts2); */
  return;
}

/*******************************************************************************
          AderivFwrts2: Compute the derivative of the CDF with respect to
                        sigma squared at x using the numerical derivative
                        ds2cum[...
*******************************************************************************/

  void AderivFwrts2()
{
  double cuminc, xpos, xposinc, smallnum;
  int    index1;
  smallnum = exp(-1.0E12);

  xpos = fabs(x);
  xpos /= sigma * twosqrt;
  if (xpos >= 20.0)
    {
     xpos = 20.0 - smallnum;
    }
 
  index1   = xpos * 10000.0;
  if (index1 >= 200000) index1 = 199999;
  cuminc   = ds2cum[index1 + 1] - ds2cum[index1];
  xposinc  = xpos - (((double)(index1))/10000.0);
  AdFwrts2 = ds2cum[index1] + ((xposinc) * cuminc);
  if (x < 0.0)
    {
     AdFwrts2 = ds2cummax - AdFwrts2;
    }
  else {AdFwrts2 += ds2cummax;}

  AdFwrts2 /= sigma2;
  ndevcdf();
  AdFwrts2 -= ndcdf / (two * sigma2);

  return;
}

/*******************************************************************************
          ds2cuminit : Compute the numerical integral of t^2 * exp[t^2] for use
                       in computing the derivative of the log-likelihood
                       with respect to sigma-squared.
*******************************************************************************/

  void ds2cuminit()
{
  double v1, v2;
  int i;

  ds2cum[0] = 0.0;

  x = 0.0001; v1 = zero;

  for (i=1; i<=200000; i++)
     {
      v2 = ( (x * x) * exp(-(x * x)) ) / pisqrt;
      ds2cum[i] = ds2cum[i-1] + ((v1 + v2)/2.0) * 0.0001;
      v1 = v2; x += 0.0001;
     }
  ds2cummax = ds2cum[200000];

printf("ds2cummax: %12.9lf\n", ds2cummax);
  return;
}

/*******************************************************************************
       alternativeA : Compute an alternative value for the derivative of
                      F(0sigma) with respect to sigma squared divided by
                      F(0,sigma) at an extreme value of x.
*******************************************************************************/

 void alternativeA()
{
  double a, b, c, d, e, xsave;

  xsave = x;
  x = 3.0;
  c = (log(x)/(sigma2 * sigma2)) + (one / (x * x * 4.0 * sigma * sigma2));
  d = one / (two * x * x);
  x = -3.0;
  ndevcdf(); derivFwrts2(); e = exp(log(fabs(dFwrts2)) - log(fabs(ndcdf)));
  a = d*e /  c;
  b = ndcdf / d;
  x = xsave;
  c = (log(fabs(x))/(sigma2 * sigma2)) + (one / (x * x * 4.0 * sigma * sigma2));
  d = one / (two * x * x);
  altA = a*(c /  d);

  return;
}

/*******************************************************************************
       numderiva : Numerically compute the derivative of the derivative
                   of the log-likelihood with respect to an alpha parameter.
*******************************************************************************/

  void numderiva()
{
 double temp11, temp12, temp13, incrementA;

 incrementA = 0.001;
 temp11 = EPest[cEvent];
 EPest[cEvent] += incrementA; derivwrta(); temp12 = dwrta;
 EPest[cEvent] -= two * incrementA; derivwrta(); temp13=dwrta;
 numdrva = (temp12 - temp13) / (two * incrementA);
 EPest[cEvent] = temp11;

 return;
}

/*******************************************************************************
       numderivb : Numerically compute the derivative of the derivative
                   of the log-likelihood with respect to a beta parameter.
*******************************************************************************/

  void numderivb()
{
 double temp11, temp12, temp13, incrementA;

 incrementA = 0.001;
 temp11 = SPest[cSPC];
 SPest[cSPC] += incrementA; derivwrtb(); temp12 = dwrtb;
 SPest[cSPC] -= two * incrementA; derivwrtb(); temp13 = dwrtb;
 numdrvb = (temp12 - temp13) / (two * incrementA);
 SPest[cSPC] = temp11;

 return;
}

/*******************************************************************************
       numderivs2 : Numerically compute the derivative of the derivative
                    of the log-likelihood with respect to sigma2.
*******************************************************************************/

  void numderivs2()
{
 double temp21, temp22, temp23, incrementB;

 incrementB = 0.0001;
 temp21 = sigma2;

 sigma2 += incrementB; sigma=sqrt(sigma2); derivwrts2(); temp22 = dwrts2;
 sigma2 -= two * incrementB; sigma=sqrt(sigma2); derivwrts2(); temp23 = dwrts2;
 numdrvs2 = (temp22 - temp23) / (two * incrementB);

 sigma2 = temp21; sigma=sqrt(sigma2);
 return;
}

/*******************************************************************************
          likelihood : Compute the value of the likelihood function
                       for given a, b, and sigma.
*******************************************************************************/

  void likelihood()
{
  /** sum1 and sum2 are computed in derivwrts2(). **/
  sigma2=sigma*sigma;

  dlikelihood=one/pow(sigma*twopisqrt,wnumuncsr);
  dlikelihood*=exp(-sigma2*sum1);
  dlikelihood+=exp(sum2);
  return;
}

/*******************************************************************************
          loglikelihood : Compute the value of the loglikelihood function
                          for given a, b, and sigma.
*******************************************************************************/

  void loglikelihood()
{
  /** sum1 and sum2 are computed in derivwrts2(). **/
  sigma2=sigma*sigma;

  dloglikelihood  = wnumuncsr*log(one/(sigma*twopisqrt));
  dloglikelihood -= sum1*sigma2;
  dloglikelihood += sum2;
/*printf("A: %12.9lf %12.9lf\n",dloglikelihood, dloglikelihood);
printf("B: %12.9lf %12.9lf\n",dloglikelihood, dloglikelihood);
printf("C: %12.9lf %12.9lf\n",dloglikelihood, dloglikelihood); */
/*
printf("wnumuncsr*log(one/(sigma*twopisqrt)): %12.9lf \n", wnumuncsr*log(one/(sigma*twopisqrt)));
printf("-sum1*sigma2:                         %12.9lf \n", -sum1*sigma2);
printf("sum2:                                 %12.9lf \n", sum2);
*/
  return;
}

/*******************************************************************************
          SigmaGrid : Find the value of sigma which minimizes dwrts2 by
                      testing a series of values for sigma.
*******************************************************************************/

 void SigmaGrid(double gridsize)
{
  int ii;
  double sigmasave, dwrts2save, wsigma, wsigma2, sigmabest, dwrts2best;

  sigmasave = sigma; dwrts2save = dwrts2;
  sigmabest = sigma; dwrts2best = dwrts2;

  sigma = sigma - (5.0 * gridsize); if (sigma < 0.0) sigma = 0.0001;

  for (ii=1; ii<=11; ii++)
    {
printf("sigma sigmagridlevel : %12.9lf %2d\n", sigma, sigmagridlevel);
     if (ii == 6) {sigma += gridsize; continue;}
     sigma2 = sigma * sigma;
     derivwrts2();
     if (fabs(dwrts2) < fabs(dwrts2best))
       {
        sigmabest = sigma; dwrts2best = dwrts2;
printf("sigmabest, dwrts2best: %12.9lf %12.9lf\n", sigmabest, dwrts2best);
       }
     sigma += gridsize;
    }

  if (sigmagridlevel == 1)
    {
     sigma = sigmabest; sigma2 = sigma * sigma;
     sigmagridlevel++;
     gridsize /= 5.0;
     SigmaGrid(gridsize);
    }

  if (sigmagridlevel == 2)
    {sigma = sigmabest; sigma2 = sigma * sigma;}

  sigmagridlevel--;
  return;
}

/*******************************************************************************
          SigmaSign : Find the value of sigma which minimizes dwrts2 by
                      testing for a change of algebraic sign.
*******************************************************************************/

 void SigmaSign(double sincrement)
{
  int sign, psign, scount;
  double prevval, sigmasave, dwrts2save, sigmabest, dwrts2best;

  sigmasave = sigma; dwrts2save = dwrts2;
  scount = 0;
  sigma -= 11.0 * sincrement; if (sigma < 0.0) sigma = 0.001;
  sigma2 = sigma * sigma;
  derivwrts2(); prevval = dwrts2;
  if (dwrts2 >= 0.0) {psign = 1;} else {psign = -1;}

  while (scount < 100)
    {
     sigma += sincrement; sigma2 = sigma * sigma;
     derivwrts2(); if (dwrts2 >= 0.0) {sign = 1;} else {sign = -1;}
     if (psign != sign)
       {
        sincrement /= -5.0;
        sigmabest = sigma; dwrts2best = dwrts2;
/* printf("dwrts2 sigma sincrement scount: %12.9lf %12.9lf %12.9lf %3d\n", dwrts2, sigma, sincrement, scount); */
       }
     psign = sign; prevval = dwrts2; scount++;
    }

  if (fabs(dwrts2best) > fabs(dwrts2save))
    {
     sigma = sigmasave; sigma2 = sigma * sigma;
    }
  else
    {
     sigma = sigmabest; sigma2 = sigma * sigma;
    }
  return;
}

/*******************************************************************************
          NewtonRaphson : Compute the parameter estimates by iterative
                          successive approximation using the Newton-Raphson
                          method (algorithm).
*******************************************************************************/

  void NewtonRaphson()
{
  double valtemp;
  int ii, cOBSsave;


/*printf("sizeof(double) sizeof(long doule): %4d %4d\n", sizeof(double), sizeof(long double));*/



/*
  sigma = one; sigma2 = one; x = -1.8;
  for (ii=1; ii<100; ii++)
    {
     ndevcdf();
     printf("x ndcdf ndcdf ndcdf: %8.3lf %18.14E  %18.14E %18.14E\n", x, ndcdf, 0.5 + erf(x/sqrt(2.0))/2.0, sndevcdfapproxC(x));
     x -= 0.1;
    }
*/

/*
 * xt = -10.0;
  for (ii=0; ii<20; ii++)
    {
     xt += 0.5;
     printf("xt ~ndcdf: %12.9lf %12.9lf\n", xt, sndevcdfapproxB(xt));
    }
*/


/*
  sigma = 1.0; sigma2 = 1.0;
  x = 0.0; xt = x;
  for (ii=0; ii<=600; ii++)
    {
     sndevcdf(); ndevcdf();
     printf("x ndcdf xsn ndcdf-xsn: %12.9lf %18.14E %18.14E %18.14E\n", x, ndcdf, xsn, ndcdf-xsn);
     x+=0.1; xt = x;
    }
*/

/*

  sigma = 1.0; sigma2 = 1.0;
  x = 0.0; 
  for (ii=0; ii<=600; ii++)
    {
     ndevcdf(); ndevpdf();
     printf(" x ndpdf/(1-ndcdf): %12.9lf %12.9lf\n", x, ndpdf/(1.0 - ndcdf));
**
     valtemp= 0.5 + erfl(x/sqrt(2.0))/2.0;
     printf ("x ndcdfA ndcdfB ndcdfA-ndcdfB : %12.9lf %12.9e %12.9e %12.9e \n", x, ndcdf, valtemp, ndcdf - valtemp);
**
     x+=0.01;
    }
*/

  numiters=0;
  ddiv=1.8;

  while (numiters<MAXITERS)
     {
      indicator=0;

      /** Compute the event parameter estimatess. **/

      cOBS=1;
      dwrtasum=0.0; dwrtacount=0.0;
      xsum = 0.0; x2sum = 0.0; xcount = 0.0;
      while (cOBS<=obscount)
          {
           cEvent=Event[cOBS];
           EPtype = EPindic[cEvent];

           switch (EPtype)
                {
                 case 0 :    /* Not-Used Case : Excluded because there is no      */
                             /* reference to this event in any file.              */
                   cOBS++;
                   break;

                 case 1 :    /* Normal Case : At least one uncensored observation */
                   
                   cOBSsave = cOBS;
                   for (ii=1; ii<=2; ii++)
                     {
                      derivwrta(); if (fabs(dwrta) < 0.000000005) break;
                      temp18 = dwrta;
                      prevEPest[cEvent]=EPest[cEvent];
                      numderiva();
                      if (fabs(numdrva) > 0.00001)
                        {
                         EPest[cEvent] -= temp18/numdrva;
                        }
                      if (isinf(EPest[cEvent]) || isnan(EPest[cEvent]))
                        {
                         EPest[cEvent] = prevEPest[cEvent];
                        }
                     }
                   derivwrta();
                   dwrtasum+=dwrta; dwrtacount+=1.0;
                   if (fabs(dwrta) < 0.0005) indicator++;
                   cOBS=cOBSnext;
                   break;

                 case 5 :     /* Missing-Data Case : Excluded because no initial   */ 
                              /* parameter estimate was provided.                  */
                   cOBS++;
                   break;

                 case 8:     /* Censored-Data-Only Case : Simple algebraic computation of event parameter */

                   cOBS++; break;
                }
          }

      /** Compute the species/cut (SPC) parameter estimates. **/
 
      inr=1;
      dwrtbsum=0.0; dwrtbcount=0.0;
      while (inr<=obscount)
          {
           cSPC=SPC[Dindex[inr]];
           SPtype=SPindic[cSPC];

           switch (SPtype)
                {
                 case 0 :    /* Not-Used Case : Excluded because there is no      */
                             /* reference to this SPC in any file.              */
                      inr++;
                      break;

                 case 1 :    /* Normal Case : This SPC is used in the data file and  */
                             /* in the file of initial parameter estimates.          */
                   for (ii=1; ii<=2; ii++)
                     {
                      derivwrtb(); if (fabs(dwrtb) < 0.00005) break;
                      temp18 = dwrtb;
                      prevSPest[cSPC]=SPest[cSPC];
                      numderivb();
                      if (fabs(numdrvb) > 0.00001)
                        {
                         SPest[cSPC] -= temp18/numdrvb;
                        }
                      if (isinf(SPest[cSPC]) || isnan(SPest[cSPC]))
                        {
                         SPest[cSPC] = prevSPest[cSPC];
                        }
                     }
                   derivwrtb();
                   dwrtbsum+=dwrtb; dwrtbcount+=1.0;
                   if (fabs(dwrtb) < 0.0005) indicator++;
                   inr=inrnext;
                   break;

                 case 5 :     /* Missing-Data Case : Excluded because no initial   */ 
                              /* parameter estimate was provided.                  */
                      inr++;
                      break;
                }
          }
 
      /** Compute the parameter estimate for sigma. **/

      for (ii=1; ii<=2; ii++)
        {
         derivwrts2(); 
         temp18 = dwrts2;
if (fabs(temp18) > 100.0) {temp18 /= 9.0;}
         prevSigma = sigma;
         numderivs2();
/*printf("sigma temp18 numdrvs2 temp18/numdrvs2: %12.9lf %12.9lf %12.9lf %12.9lf\n\n", sigma, temp18, numdrvs2, temp18/numdrvs2); */
         if (fabs(numdrvs2) > 0.000001)
           {
            sigma2 -= temp18/numdrvs2;
            sigma   = sqrt(sigma2);
           }
         if (isinf(sigma2) || isnan(sigma2))
           {
            sigma = prevSigma;
            sigma2 = sigma * sigma;
           }
        }

/* if (sigma < 0.420 || sigma > 0.459) {sigma=0.447; sigma2=sigma*sigma;} */

/*   sigma = 0.447; sigma2 = sigma * sigma; SigmaSign(sigma / 100.0); */

/*   if (fabs(dwrts2) < 60000.0) {sigmagridlevel = 1; SigmaGrid(sigma / 300.0);}
     else printf("dwrts2: %12.9lf\n", dwrts2); */


     /** Print intermediate values. **/

      numiters++;
      derivwrts2(); loglikelihood();
      if (dloglikelihood > maxloglikelihood)
          {
           maxloglikelihood=dloglikelihood;
           maxlogiteration=numiters;
           sigmaBest=sigma;
           for (i=0; i<SIZEEVENTARRAY; i++)
               {
                if (EPindic[i]==1) EPBest[i]=EPest[i];
                if (EPindic[i]==8) EPBest[i]=EPest[i];
               }
           for (i=0; i<SIZESPCARRAY; i++)
               {
                if (SPindic[i]==1) SPBest[i]=SPest[i];
               }
          }

      printf("[%5d] %18.9lf | %10.8lf | %12.9lf // %14.9lf -- %14.9lf [%5d]\n",
             numiters, dwrts2, sigma, dloglikelihood, dwrtasum/dwrtacount, dwrtbsum/dwrtbcount, extremecount);

      /** Check for convergence. **/

      llsave[numiters] = dloglikelihood;
      if (indicator==eventcount+spccount+1)
          {
           if (numiters > 25 && fabs(llsave[numiters-20] - maxloglikelihood) < 0.000000001) break;
          }
      if ((dwrtasum/dwrtacount < 0.0005 && dwrtbsum/dwrtbcount < 0.0005 && dwrts2 < 0.0000001))
          {
           if (numiters > 25 && fabs(llsave[numiters-20] - maxloglikelihood) < 0.0000000001) break;
          }

     } /** End While Loop for iterations. **/


  return;
}
