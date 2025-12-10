# ML Project: Single-cell Classification in Nematodes

This project focuses on the analysis and classification of nematode single-cell data
based on gene expression profiles. The goal is to explore unsupervised structure in
the data and build machine learning models for cell-type classification.

## Project structure

1. **Background / Introduction**  
   Short overview of the biological context and the dataset.

2. **Data exploration**  
   Basic preprocessing, quality checks and exploratory analysis of single-cell
   gene expression data.

3. **Clustering**  
   Dimensionality reduction and clustering (e.g. Leiden) to identify putative
   cell populations in an unsupervised way.

4. **Classical machine learning models**  
   Training and evaluation of several scikit-learn models based on the
   extracted features / embeddings:
   - Support Vector Machines (SVM)
   - Random Forest
   - K-Nearest Neighbours (KNeighborsClassifier)
   - Multinomial Logistic Regression
   - MLPClassifier

5. **Neural networks \& Transformers**  
   Experiments with deeper neural network architectures, including
   multi-layer perceptrons and simple multi-head Transformer models
   applied to the gene expression representations.

6. **Possible extensions**  
   Ideas for further work, e.g. better feature engineering, improved
   architectures, benchmarking on additional datasets.

## Methods & tools

- **Language & libraries:** Python, scikit-learn, PyTorch
- **Tasks:** clustering, classification, model comparison
- **Environment:** originally developed in a Kaggle notebook and adapted
  to this repository for demonstration and reproducibility.
