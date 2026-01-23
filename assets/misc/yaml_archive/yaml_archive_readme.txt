# Here lie previous YAML version for weighting gene predictions from Ragnarok.
# For the best results, please ensure you are using the `plant_microexon_filter.yaml` file.

# 1. [Deprecated] plant.yaml is the default plant yaml from mikado.
# 2. [Deprecated] plant_min_exon_length.yaml applies a 3-bp hard filter to remove microexon genes.
# 3. [Current] plant_microexon_filter.yaml applies a 1.5x negative multiplier per basepair to microexons. The most severe penalty are for exons of length 1, and the least severe are for exons of length 15. The maximum penalty for a model with a microexon length of 1 is -22.5.
