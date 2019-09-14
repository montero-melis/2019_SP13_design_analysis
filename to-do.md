TO DO -- SP13 design analysis
============================

bfda_sp13
---------

- State B=10000 in the intro!
- Make up our mind about Chen et al's guidelines vs formula - leave only corresponding simulations

reanalysis_original
-------------------

- Section on effect size - say in what sense this is small
- Remove 2x4 ANOVA analysis from re-analysis script (to make it less cumbersome)?


Analysis pipeline
-----------------

- Generate data with a small effect (use set.seed())
- Compute "standard" BF in the way we will for the actual study and make sure it is inconclusive with the first 60 participants
- Then generate data from 12 additional participants, compute BF again
- Make sure that the evidence threshold (BF>6) is reached before N=96
- In the end compute both the standard and the replication BF and report.
