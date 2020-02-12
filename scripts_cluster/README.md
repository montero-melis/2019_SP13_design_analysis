README
======

This folder was created to run the scripts on a computer cluster so that the
simulations could be sped up. It works independently of other files in the
repo and that is why there is a bit of duplication. So for instance the
script "generate_data_fnc.R" is a duplicate from that in the root directory
of the repo. And similarly for the script that runs the simulations: much of
it is a duplication of the knitr document "Appendix_C_power-analysis_sp13.Rmd".

This could lead to some things not working properly if that's not taken into
account (e.g., if you run the script "power-simulations_sp13_par.R") from the
Rproject root directory it will assume the root path and will therefore lead
to error messages. Keep that in mind!
