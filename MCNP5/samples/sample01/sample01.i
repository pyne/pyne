DAG-MCNP5 Sample Problem 01 - water sphere
C Data cards
C ********************
C --------------------
C Materials
C --------------------
C Material 1: Water
M1 1001 2 8016 1
C --------------------
C Source
C --------------------
C isotropic 14 MeV point source at the origin
SDEF
C --------------------
C Job Control
C --------------------
C run in neutron only mode
mode n
C run for 1000 histories
nps 1000
