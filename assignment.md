### 1. Assignment - Position Weight Matrices

- how many candidates are there for the different start codon variants?
- find a threshold for the logarithmic score such that 50% of the valid candidates are detected
    - parameters: lenght: 30, background_dist: uniform, pseudo_count: 1
    - how many false positives?
- divide data into training (first 400 sequences) and test set (rest)
- use training set to determine threshold needed to classify 50% of the valid ones right
