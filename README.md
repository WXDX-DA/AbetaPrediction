# Integrated Algorithm Combining Plasma Biomarkers and Cognitive Assessments Accurately Predicts Brain β-Amyloid Pathology
This repository contains the source code for the research paper "Integrated algorithm combining plasma biomarkers and cognitive assessments accurately predicts brain β-amyloid pathology" published in Communications Medicine. The code provided here is open-source and can be used to reproduce the results and analyses described in the paper. 

[![DOI](https://zenodo.org/badge/626207371.svg)](https://zenodo.org/badge/latestdoi/626207371)

## Overview
The research paper presents a novel algorithm that combines plasma biomarkers and cognitive assessments to predict brain β-amyloid pathology. The algorithm is designed to improve the early detection and monitoring of Alzheimer's disease, which is characterized by the accumulation of β-amyloid plaques in the brain. This repository includes the code for preprocessing data, training the machine learning models, and evaluating the performance of the algorithm.

## Requirements
To install and run the code, you will need the following:

* R version 4.20 or higher
* tidyverse
* foreign
* kableExtra
* pROC
* caret
* rpart
* rpart.plot
* gtsummary
* scales
    
You can install the required libraries using the following command in R:

````
install.packages(c("tidyverse", "foreign", "kableExtra", "pROC", "caret", "rpart", "rpart.plot", "gtsummary", "scales"))
````

## Repository Structure
The repository is organized as follows:

  * Demographic_DementiaPrediction.r: R script for demographic analysis and dementia prediction.
  * Tree_performances-zscore.r: R script for evaluating decision tree performance using z-score transformed data.
  * ADNI_validation.r: R script for validating the algorithm using the ADNI dataset.
  * Single_var_prediction.r: R script for single-variable prediction.
  * PETwithHistory.r: R script for analyzing PET scans combined with patient history.

## Usage
1. Clone the repository:
```
git clone https://github.com/WXDX-DA/AbetaPrediction.git
cd AbetaPrediction
```
2. Install the required libraries in R:

```
install.packages(c("tidyverse", "foreign", "kableExtra", "pROC", "caret", "rpart", "rpart.plot", "gtsummary", "scales"))
```

3. Open RStudio or your preferred R environment, and run the provided R scripts:

* Demographic analysis and dementia prediction:
```
source("Demographic_DementiaPrediction.r")
```
* Evaluating decision tree performance using z-score:
```
source("Tree_performances-zscore.r")
```
* Validating the algorithm using the ADNI dataset:
```
source("ADNI_validation.r")
```
* Single-variable prediction:
```
source("Single_var_prediction.r")
```
* Analyzing PET scans combined with patient history:
```
source("PETwithHistory.r")
```

## Citation
If you use this code or the findings from the research paper, please cite the original publication.


## License
This project is licensed under the GNU GPLv3 License - see the [LICENSE](https://github.com/WXDX-DA/AbetaPrediction/blob/main/LICENSE) file for details.

## Contact
For questions or support, please contact the corresponding author of the paper or raise an issue on this GitHub repository.
