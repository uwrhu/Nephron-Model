# Nephron-Model
This project is building sex-specific computational models for the epithelial transport along rat nephron.

If Python3 is default, run ‘python simulation.py --sex Male/Female --diabete Y/N --inhibition NHE3/NKCC2/NCC/ENaC/SNB --percentage 0.5/0.7/1.0'.
If not, replace 'python' with 'python3'.

Options:

Sex: Male or Female, required.

diabete: Y or N, required.

inhibition: NHE3, NKCC2, NCC, ENaC, SNB(sequential nehpron blockade), optional. (No inhibition means normal nephron).

All the output files’ names are in following structure: ‘sex+segment’_’concentration or flow’_of_’solute’_in_’compartment’.txt.

Here is an example: femaleccd_con_of_Cl_in_Bath.txt. It contains interstitial concentration of Chloride along cortical collecting duct in female rat.

Another example: malept_flow_of_Na_in_Lumen.txt. It contains luminal flow of Sodium along proximal convolute tubule in male rat.

These results are scaled per nephron.

The unit of concentration from outputs is mmol/L (mM).

The unit of volume is nl/min.

The unit of flow is pmol/min.
