# Explanation of PSI Adjustment for Exon Skipped Reads

## Addressing Batch Effects in Splicing Reporter Analysis

In our analysis, we encounter a phenomenon known as "PSI depression." This occurs when sequences that are expected to be consistently included in the analysis are not, due to factors such as PCR cycles and the presence of chimeric reads. To address this, we identify ~100 "anchor" sequences that should ideally have a Percent Spliced In (PSI) value of 1, indicating perfect inclusion. These sequences are selected by calculating the geometric mean of PSI values for each element across all samples and choosing the top 100 sequences with the highest geometric mean.

## Modeling the Barcode Swapping Rate

To understand the impact of barcode swapping, we define the swapping rate as $x$, and the total read count as $T$. The observed number of exon-skipped reads for each barcode sequence is denoted as $\hat{E}$. This observed value is calculated by taking the true exon-skipped reads, subtracting the reads that swapped to a different barcode, and adding the reads that swapped to the barcode:

$$
\hat{E} = E - xE\frac{T-E}{T} + x(T-E)\frac{E+I}{T}.
$$

Given the large number of sequences, where each contributes minimally to the total reads, we assume $T \gg E$. This assumption simplifies the equation to:

$$
\begin{align*}
\hat{E} &= E - xE\frac{T-E}{T} + x(T-E)\frac{E+I}{T} \\
\hat{E} &= E - xE + x(E+I) \\
\hat{E} &= E + xI 
\end{align*}
$$

The observed PSI, denoted as $\hat{\text{PSI}}$, is calculated as $\frac{I}{I + \hat{E}}$. The true PSI is given by $\text{PSI} = \frac{I}{I+E}$. In the ideal scenario where $\text{PSI} = 1$, the exon-skipped reads $E$ are zero. By substituting into the PSI equation, we can solve for the swapping rate $x$:

$$
\begin{align*}
\hat{\text{PSI}} &= \frac{I}{I + E + xI} \\
x &= \frac{1-\hat{\text{PSI}}}{\hat{\text{PSI}}}
\end{align*}
$$

To correct for this, we calculate the swapping rate $x$ for each of the ~100 anchor sequences and then take the average. Using this average swapping rate, we can adjust the actual exon-skipped reads $E$ for each sequence and subsequently compute the true PSI values.
