# PreCDA
## Introduction
Circular RNAs (circRNAs) are a novel class of endogenous noncoding RNAs. Emerging evidence has shown that circRNAs can be novel biomarkers or therapeutic targets for many diseases. Therefore, identifying potential disease-related circRNAs is helpful in improving the efficiency of finding therapeutic targets for diseases. Here, a computational model (PreCDA) is proposed to predict potential circRNA-disease associations.
## Usage
IDE：IntelliJ IDEA  
Development Language: Java, scala  
Note: please set the appropriate path before running

Step 1. run filehandler.java to create files taht contains circRNA-disease associations from different circRNA databases  
Step 2. run ExpressProfile.java to calculate circRNA expression similarity  
Step 3. run circRNAsim.java to calculate circRNA functional similarity  
Step 4. run UniNet.java to build the associations between circRNAs  
Step 5. run PersonalRank.java to rank candidate circRNAs  
Step 6. run TestSet.java to output the AUC value of predicting candidate disease-related circRNAs 
