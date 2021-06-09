## Procedure for explicit local (local time stepping - LTS/lts) w/ IBM

1) Get initial solution
   This can be obtained through either a RANS simulation or performing a single iteration w/o LTS
2) prep_lts_ibm.py w/ the t.cgns file obtained from (1) and tc.cgns needed to perform (1)
3) compute.py
