/* Heart failure survival analysis - SAS version
   Data: Heart Failure Clinical Records (UCI), 299 patients
   Purpose: repeat the R analysis in SAS and check that the results match */

/* Put the path to the CSV file here, in one place */
%let csvpath = /home/u64544172/heart_failure_clinical_records_dataset.csv;

/* ---- Read the data ---- */
proc import datafile="&csvpath"
    out=work.hf
    dbms=csv
    replace;
    getnames=yes;
    guessingrows=max;
run;

/* Check: 299 observations, 13 variables */
proc contents data=work.hf;
run;

/* DEATH_EVENT: 1 = died, 0 = alive at the end (censored) */
proc freq data=work.hf;
    tables DEATH_EVENT;
run;

/* ---- Kaplan-Meier: survival for all patients ---- */
proc lifetest data=work.hf method=km timelist=50 100 150 200 250 reduceout;
    time time*DEATH_EVENT(0);
run;

/* ---- Two groups by ejection fraction ---- */
data work.hf2;
    set work.hf;
    length ef_group $14;
    if ejection_fraction < 40 then ef_group = "EF below 40";
    else ef_group = "EF 40 or above";
run;

proc freq data=work.hf2;
    tables ef_group*DEATH_EVENT;
run;

/* ---- Log-rank test: do the two groups survive differently? ---- */
/* notable hides the long list of every patient time */
proc lifetest data=work.hf2 method=km notable;
    time time*DEATH_EVENT(0);
    strata ef_group / test=logrank;
run;

/* ---- Cox regression ---- */
/* ties=efron gives the same results as R's coxph() */
proc phreg data=work.hf2;
    model time*DEATH_EVENT(0) = age ejection_fraction serum_creatinine
          serum_sodium anaemia high_blood_pressure diabetes sex smoking
          / risklimits ties=efron;
run;