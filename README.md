# Guelph-Family-Health-Study-published-
Background:
Taste is a fundamental determinant of food selection, and inter-individual variations in taste perception may be important risk factors for poor eating habits and obesity. Characterizing differences in taste perception and their influences on dietary intake may lead to an improved understanding of obesity risk and a potential to develop personalized nutrition 
recommendations. 

Study process:
This study explored associations between 93 single nucleotide polymorphisms (SNPs) in sweet, fat, bitter, salt, sour, and umami taste receptors and psychophysical measures of taste. Forty-four families from the Guelph Family Health Study participated, including 60 children and 65 adults. Saliva was collected for genetic analysis and parents completed a three-day food record for their children. Parents underwent a test for suprathreshold sensitivity (ST) and taste preference (PR) for sweet, fat, salt, umami, and sour as well as a phenylthiocarbamide (PTC) taste status test. Children underwent PR tests and a PTC taste status test.

Result and impact:
Analysis of SNPs and psychophysical measures of taste yielded 23 significant associations in parents and 11 in children. After adjusting for multiple hypothesis testing, the rs713598 in the TAS2R38 bitter taste receptor gene and rs236514 in the KCNJ2 sour taste-associated gene remained significantly associated with PTC ST and sour PR in parents, respectively. In children, rs173135 in KCNJ2 and rs4790522 in the TRPV1 salt taste-associated gene remained significantly associated with sour and salt taste PRs, respectively. A multiple trait analysis of PR and nutrient composition of diet in the children revealed that rs9701796 in the TAS1R2 sweet taste receptor gene was associated with both sweet PR and percent energy from added sugar in the diet. These findings provide evidence that for bitter, sour, salt, and sweet taste, certain genetic variants are associated with taste function and may be implicated in eating patterns. (Support was provided by the Ontario Ministry of Agriculture, Food, and Rural Affairs).

Summary:
The statistical analysis of this Guelph family health study consists 2 parts: 1: taste preference analysis 2. diet composition analysis
Variables of interest: taste receptors, measure of taste; diet variables; demographics variables (as confounding variables); 

Breakthrough:1: family data are correlated as clusters; pedigree relationship is given as ID of individual ID and parents IDs
             2: Confounding variables controlled
             
Part 1: taste preference and sensitivity threshold

Analysis 1: inverse mapping: SNP ~ age+sex+BMI+taste variable  using GEE
             mapping: taste ~ SNP+age+sex+BMI    using GEE  
         2: SNP ~ taste using GQLS (can only do one to one)
            SNP ~ sex, age, BMI, ethics using GQLS
            
         we found that sex, age, BMI, ethics are confounding variables, we need to incorporate them in GEE model
         
Part 2: diet composition and taste

Analysis 1: diet ~ age+sex+BMI
         2: SNP ~ taste+diet
      

            
          
            
