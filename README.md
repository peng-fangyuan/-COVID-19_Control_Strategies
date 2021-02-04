# COVID-19_Control_Strategies
Data and Code for Control Strategies against COVID-19 in China: Significance of Effective Testing in the Long Run


## Code

### Inference Code
We obtain essential parameters by mainly using inference.m. 
SEIR.m and initial.m will be used in inference.m.

### Counterfactual Analysis Code

After obtaining inferred parameters, we conduct scenarious analysis by using  counter_facutal_t.m

## Data

1. Intercity Mobility Data
2020: Mig_post.mat
2019: Mig2019_post.mat

2. Within-city Mobility Data
2020: in2020_post.mat
2019: in2019_post.mat

3. Confirmed case data 
incidence.mat 

4. Population data 
pop.mat

5. City name data 
city.mat

6. Hubei city 
hu.mat is to check which city is in Hubei. 

Note: Observed data is from Jan 1 to Mar 15. For mobility data, the data after Mar 15 is the mean value for each city of the week before Mar 15, which use to conduct long term counterfacutal analysis.   
