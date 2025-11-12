12/11/2025 Update:
  Change the method from IBM simulation to Markov Chain with Reward.
  Fix the graph problem so it is more closed to the analytical results of mean/variance/skewness of lifespan/LRO.
  Questions remained to be solved:
  1. How to combine the results of different species together and to interpret that result?



28/10/25 Update:
  Deleted the IPM part as it's unnecessary. Updated the Simulate IBM-GT.R file: 1. add no. of survived chicks/female to the model to make it more realistic 2. add mother_ID to the sim.data to make further calculation of LRO easier (probably)
  Questions remained to be solved:
  1. How to do the LRO and reproductive senescence graph -- also use sim.data. Finish by 29/10/25
  2. Add counterfactual model into the codes! This might require another R.script for counterfactual models.




26/10/25 Update:
  Upload the code for the test run. 
  Test dataset: Great Tit from Jones et al., 2014.
  Test Method: Use IBM to simulate populations and randomly select individuals. While use the survival and fecundity derived from the dataset to build IPM. Then compare the difference in lifespan distribution and fecundity between IBM and IPM.
  Questions remained to be solved: 
    1. The graph compares IBM model VS IPM model (both contains senescence), instead of conterfactual model VS +senescence model.
    2. The fecundity graph just describes how fecundity changes with age instead of how LRO changes with age/lifespan.
    3. How to interpret the two final graphs? There is difference between different models, but how could we interpret the difference? Can we run a stats test for the difference?
