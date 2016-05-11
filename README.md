QTLSimulator
-----------

QTL Simulator simulates QTLs with c number of celltypes with for every QTL and every sample:

```
simulatedExpression = 0
for number of celltypes do:
    simulatedExpression += (B1 * cellcounts) + (B2 * genotype);
done
```

where

    B1 = cellcount coefficient, e.g. if B1 > 0 the expression is correlated to cellcount  
    B2 = genotype coefficient (the QTL effecT). If B2 > 0 genotype has positive effect on expression, B2 < 0 genotype has negative effect on expression
    genotype = Random number 0, 1, or 2
    cellcounts = Normal distribution with mean a percentage given by the user (with commandline option -c)

This is done for -n QTLs with different genotype coefficient groups. The output is 4 files, a simulated expression, genotype, and cellcount file, and an info file with the B1 and B2 per QTL


Run
---
java -jar QtlSimulator.jar -c 70,20,7,3 -n 100, -o testOutFolder/ -s 100 -e 4

c: Cellcount percentages  
n: Number of QTLs to write  
o: Outfolder  
s: Samplesize  
e: Noise  

