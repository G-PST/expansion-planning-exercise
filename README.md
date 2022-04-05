# expansion-planning-exercise

This repository demonstrates a coordinated workflow using several different tools. The modeling steps are as follows:

1. Start with [RTS-GMLC](reference-matpower/RTS_GMLC/)
2. PyPSA
   1. Adjust system requirement
      1. Remove conventional generation (remove coal) / Set RE target
      2. Increase load (reflecting standard load growth, electrification, ...)
   2. Optimal expansion planning (PyPSA)
      1. Assumptions: only add/adjust conductors on existing paths, capital-costs, ...
      2. One scenario results in a new generation portfolio
      3. One scenario results in a new generation and transmission porfolio
3. Optimal integration planning (Powsybl, Pandapower)
   1. establish assumptions for representing details of expansions identified in step 3. (e.g. if 3.3 identifies new transmission capacity along a path, should we optimize the design of the entire path or just the new capacity?)
   2. enforce DC security constraints
   3. results in a N-1 secure, AC feasible (hopefully) network
4. System performance (no new resource)
   1. Load flow
   2. UC/ED analysis
   3. Production Cost Modeling
