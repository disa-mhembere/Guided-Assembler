Guided Assembler Game plan
=========================

Data Structures:
----------------
1. Get genome for 2 `simple` and `similar` organisms (ask about prokaryotic organisms) 
2. Take one, split up into reads, add errors and use as synthetic data 
3. Code:
  - *BWT & accompanying updater (Dynamic BWT) --> DM*
  - *Suffix Array & updater (Dynamic SA) --> DM*
  - Dynamic Tally Array - we don`t how to do yet --> DM (SL)
4. Aligner:
  - Simple aligner (See lecture 12, slide 7ish) --> SL

Timeline: All DS-related code check-in Fiday November 22nd. 

Assembly algorithm:
------------------
Timeline: Start: Fri November 22nd - End: Th November 28th code check-in
Generally:
- Exact match for reads via aligner
- *Fill in gaps/update reference*

Testing:
-------
- Synthtic eval i.e compare assemled genome with reference (single genome)
- Use both genomes as references against the other to attempt reconstruct. (dual genome). This is 2 experiments.

Write up & Presentation:
------------------------
TODO: 


