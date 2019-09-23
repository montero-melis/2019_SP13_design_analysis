TO DO -- SP13 design analysis
============================

reanalysis_original
-------------------

- mean effects using "emmeans" package (Phillip A's suggestion)?


Analysis pipeline
-----------------

- Generate data with a small effect (use set.seed())
- Compute "standard" BF in the way we will for the actual study and make sure it is inconclusive with the first 60 participants
- Then generate data from 12 additional participants, compute BF again
- Make sure that the evidence threshold (BF>6) is reached before N=96
- In the end compute both the standard and the replication BF and report.
