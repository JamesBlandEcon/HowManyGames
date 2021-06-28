This folder contains:

 - HowManyGames.csv, and HowManyGames.dta: Data from the experiment
 - MATLABdata.mat: Data stored in Matlab format
 - HowManyGames.do: Stata .do file that produces the regression tables reported in the paper. These files are exported into the folder "outputs"
 - Plots.m: Matlab .m file that produces the plots in Section 3

Structural model files (Matlab)
 - MixmodelBayes2.m: Simulates and stores the posterior distribution of the model heterogeneous agents model (i.e. run this first)
 - Analyzebayes.m: Postsimulation analysis of the heterogeneous agents structural model (i.e. run this second)
 - QRE_empiricalMH_split.m: main file for estimating the QRE models
 - LearningModelPredictionsCheck.m: group mean predictions from the learning model (Figure 9)
 - All other files are functions that are called by one of the above scripts.
You may need to change some of the file paths to ones that exist on your machine.

Data in HowManyGames.dta:
  obs:         2,560                          
 vars:            66                          28 Nov 2017 12:35
 size:       360,960                          
---------------------------------------------------------------------------------------------------------------------------------------
              storage   display    value
variable name   type    format     label      variable label
---------------------------------------------------------------------------------------------------------------------------------------
Date            long    %12.0g                 Date
a               byte    %8.0g                 Payoff in Game 1 if volunteer
c               int     %8.0g                 Payoff in Game 2 if volunteer
x               int     %8.0g                 Payoff in Game 1 if free ride
y               int     %8.0g                 Payoff in Game 2 if free ride
Treatment       byte    %8.0g                  Treatment
d1              byte    %8.0g                 Points added to Game 1 relative to Treatment 1
d2              byte    %8.0g                 Points added to Game 2 relative to Treatment 1
chi1_1          float   %9.0g                 multiplier in Lottery Task 1.1
chi1_2          float   %9.0g                 multiplier in Lottery Task 1.2
chi1_3          float   %9.0g                 multiplier in Lottery Task 1.3
chi1_4          float   %9.0g                 multiplier in Lottery Task 1.4
chi2_1          float   %9.0g                 multiplier in Lottery Task 2.1
chi2_2          float   %9.0g                 multiplier in Lottery Task 2.2
chi2_3          byte    %8.0g                 multiplier in Lottery Task 2.3
chi2_4          byte    %8.0g                 multiplier in Lottery Task 2.4
ugroupid        long    %12.0g                Unique identifier for groups
uid             long    %12.0g                Unique identifier for subjects
upartnerid      long    %12.0g                Partner's uid
Session         byte    %8.0g                  Session
Period          byte    %8.0g                  Period
Subject         byte    %8.0g                  Subject
Group           byte    %8.0g                  Group
ActionAB        byte    %8.0g                 Action in Game 1, 1==B
ActionCD        byte    %8.0g                 Action in Game 2, 1==D
NumAC           byte    %8.0g                 Number of times action AC played in previous period
NumAD           byte    %8.0g                 Number of times action AD played in previous period
NumBC           byte    %8.0g                 Number of times action BC played in previous period
NumBD           byte    %8.0g                 Number of times action BD played in previous period
CountAC         byte    %8.0g                 Number of times action AC played in group in all past periods
CountAD         byte    %8.0g                 Number of times action AD played in group in all past periods
CountBC         byte    %8.0g                 Number of times action BC played in group in all past periods
CountBD         byte    %8.0g                 Number of times action BD played in group in all past periods
Partner         byte    %8.0g                  Partner
partnerAB       byte    %8.0g                 Partner's action in Game 1, 1==B
partnerCD       byte    %8.0g                 Partner's action in Game 2, 1==D
NumCorrectNGQ   str13   %13s                   NumCorrectNGQ
correctNGQ_1    byte    %8.0g                 Answered game question 1 correctly
correctNGQ_2    byte    %8.0g                 Answered game question 2 correctly
correctNGQ_3    byte    %8.0g                 Answered game question 3 correctly
correctNGQ_4    byte    %8.0g                 Answered game question 4 correctly
correctNGQ_5    byte    %8.0g                 Answered game question 5 correctly
correctNGQ_6    byte    %8.0g                 Answered game question 6 correctly
correctNGQ_7    byte    %8.0g                 Answered game question 7 correctly
invest1_1       byte    %8.0g                 Tokens invested in Lottery Task 1.1
invest1_2       byte    %8.0g                 Tokens invested in Lottery Task 1.2
invest1_3       byte    %8.0g                 Tokens invested in Lottery Task 1.3
invest1_4       byte    %8.0g                 Tokens invested in Lottery Task 1.4
invest2_1_1     byte    %8.0g                 Red tokens invested in Lottery Task 2.1
invest2_1_2     byte    %8.0g                 Red tokens invested in Lottery Task 2.2
invest2_1_3     byte    %8.0g                 Red tokens invested in Lottery Task 2.3
invest2_1_4     byte    %8.0g                 Red tokens invested in Lottery Task 2.4
invest2_2_1     byte    %8.0g                 Blue tokens invested in Lottery Task 2.1
invest2_2_2     byte    %8.0g                 Blue tokens invested in Lottery Task 2.2
invest2_2_3     byte    %8.0g                 Blue tokens invested in Lottery Task 2.3
invest2_2_4     byte    %8.0g                 Blue tokens invested in Lottery Task 2.4
NumCorrectLT1Q  str13   %13s                   NumCorrectLT1Q
correctLT1Q_1   byte    %8.0g                 Answered lottery task part 1 question 1 correctly
correctLT1Q_2   byte    %8.0g                 Answered lottery task part 1 question 2 correctly
correctLT1Q_3   byte    %8.0g                 Answered lottery task part 1 question 3 correctly
NumCorrectLT2Q  str13   %13s                   NumCorrectLT2Q
correctLT2Q_1   byte    %8.0g                 Answered lottery task part 2 question 1 correctly
correctLT2Q_2   byte    %8.0g                 Answered lottery task part 2 question 2 correctly             
---------------------------------------------------------------------------------------------------------------------------------------

