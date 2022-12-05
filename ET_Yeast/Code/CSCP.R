#devtools::install_github("rdinnager/slimr")
library(tidyverse)
library(cowplot)
library(slimr)

slim_script(
  slim_block(
    initialize(),
    {
      setSeed(12345);
      initializeMutationRate(1e-7);
      initializeMutationType("m1", 0.5, "f", 0.0);
      initializeGenomicElementType("g1", m1, 1.0);
      scriptForQTLs = "rexp(1);";
      initializeMutationType("m2", 0.5, "s", scriptForQTLs);
      initializeGenomicElementType("g2", m2, 1.0);
      m2.convertToSubstitution = F;
      m2.mutationStackGroup = -1;
      m2.mutationStackPolicy = "l";
      m1.color = "green";
      m2.color = "red";
      
      defineConstant("C", 100);
      defineConstant("W", 5000);
      pos = 0;
      q = NULL;
      
      for (i in 1:C)
        initializeGenomicElement(g1, pos, pos + W-1);
        pos = pos + W;
        initializeGenomicElement(g2, pos, pos);
        q = c(q, pos);
        pos = pos + 1;
        initializeGenomicElement(g1, pos, pos + W-1);
        pos = pos + W;
        
      defineConstant("Q", q);
      initializeRecombinationRate(1e-8);
    }),
  slim_block(1, early(),
             {
               sim.addSubpop("p1", 10000);
               writeFile(paste0("/Users/ezra/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/genome_CSCP",".csv"), paste0("Generation,Position,Frequency,Effect, origin"));
               
               writeFile(paste0("/Users/ezra/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/MeanPhenotypes_CSCP",".csv"), paste0("Generation,Phenotype"));
             }),
  
  #slim_block(mutationEffect(m2),
   #          {
  #             return 1.0;
   #          }),
     
  
  slim_block(1,..,late(),
             {
               inds = sim.subpopulations.individuals;
               phenotype = inds.sumOfMutationsOfType(m2);
               optimum = 10.0;
               inds.fitnessScaling = dnorm(optimum - phenotype, 0.0, 1.0);
               inds.tagF = phenotype;
             }),
  slim_block(1,1000,late(),
             {
               if (sim.cycle == 1)
                 cat("mean phenotype:\n");
               meanPhenotype = mean(p1.individuals.tagF);
               cat(format("%.2f", meanPhenotype));
               
               cat("\n\n-------------------------------\n");
               cat("QTLs at cycle " + sim.cycle + ":\n\n");
               
               qtls = sim.mutationsOfType(m2);
               freq = sim.mutationFrequencies(NULL, qtls);
               selCoef = qtls.selectionCoeff;
               posit = qtls.position;
               origin = qtls.originTick;
               indices = order(freq, F);
               
               for (i in indices)
                 writeFile(paste0("/Users/ezra/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/genome_CSCP",".csv"), paste0(sim.cycle,",", posit[i], ",", format("%.6f",freq[i]), ",", format("%.6f",selCoef[i]), ",", origin[i]), append = T);
                 meanPhenotype = mean(p1.individuals.tagF);
                 writeFile(paste0("/Users/ezra/OneDrive - University of Missouri/Academics Mizzou (etb68@umsystem.edu)/Biologycal_Sciences/Research_Projects/Yeast_GitHub/evogen-sims/ET_Yeast/output.dir/MeanPhenotypes_CSCP",".csv"), paste0(sim.cycle,",",format("%.6f", meanPhenotype)), append=T);
             }),
  slim_block(1000,late(),
             {
               cat("\n\n-------------------------------\n");
               cat("QTLs at generation " + sim.cycle + ":\n\n");
               
               qtls = sim.mutationsOfType(m2);
               freq = sim.mutationFrequencies(NULL, qtls);
               selCoef = qtls.selectionCoeff;
               posit = qtls.position;
               origin = qtls.originTick;
               indices = order(freq, F);
               
               for (i in indices)
                 cat("   " + posit[i] + ": selCoef = " + selCoef[i] + ", freq == " + freq[i] + ", origin == " + origin[i] + "\n");
               sim.simulationFinished();
             })
) -> def_slim1

def_slim1

slim_script_render(def_slim1)


