# Explanation of PSI Adjustment for Exon Skipped Reads - addressing batch effects from barcode swapping


In our analysis, we encounter a phenomenon known as "PSI depression." This occurs when sequences that are expected to be consistently included in the analysis are not, due to factors such as PCR cycles and the presence of chimeric reads. 

## Version 1: Modeling the Barcode Swapping Rate using PSI

To address this, we identify ~100 "anchor" sequences that should ideally have a Percent Spliced In (PSI) value of 1, indicating perfect inclusion. These sequences are selected by calculating the geometric mean of PSI values for each element across all samples and choosing the top 100 sequences with the highest geometric mean.

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

## Version 2: Modeling the barcode swapping rate from "chimeric reads" percentage from FASTQ alignment on included reads

Another way we can do it is to model the barcode swapping rate from the percentage of "chimeric reads" from the FASTQ alignment on the included reads. This is a metric that we calculate already when we do the mapping - where for all reads that are included, we try to see if they map to the sequence that's defined by the barcode at the end. If it's not a match, we call it as a "chimeric read".

We can use this metric to estimate the barcode swapping rate by assuming that the number of chimeric reads is proportional to the number of reads that swapped to a different barcode. 

Say the chimeric rate is $C$. For each exon-skipped element, the number of chimeric read is $y$. The observed total number of skipped reads is $\hat{E} = E + y$. The true number of skipped reads is $E$. So the chimeric rate $C$ is $\frac{y}{E + y}$. So we have the following equations:

$$
\begin{align*}
\hat{E} &= E + y \\
C &= \frac{y}{E + y}
\end{align*}
$$

Solving for $E$ and $y$, we get:

$$
\begin{align*}
y = \hat{E}C \\
E = \hat{E} - y = \hat{E} - \hat{E}C = \hat{E}(1-C)
\end{align*}
$$

## Version 3: Also including the included reads

I feel like version 2 is not good enough (or not as comparable, beause it matters how many included reads there are). So we can also include the included reads in the calculation.

Say there are $I$ included reads and $E$ exon-skipped reads. The chimeric rate is $C$. Then the total reads is $T = I + E$. The observed number of exon-skipped reads is $\hat{E} = E + C(I + E)$. 

Solving for $E$, we get:

$$
\begin{align*}
\hat{E} &= E + C(I + E) \\
E &= \frac{\hat{E} - CI}{1+C}
\end{align*}
$$

In the case where there are multiple E elements, then we can do the same calculation for each element and take the weighted values. E.g. if there are $n$ elements, $E_1, E_2, \ldots, E_n$, then the observed number of exon-skipped reads is $\hat{E} = E_1 + E_2 + \ldots + E_n$. The chimeric rate is $C$. Then the total reads is $T = I + E_1 + E_2 + \ldots + E_n$. The observed number of exon-skipped reads is:

$$
\hat{E} = (E_1 + E_2 + \ldots + E_n) + C(I + E_1 + E_2 + \ldots + E_n)
$$

Solving for $E_1, E_2, \ldots, E_n$, we get:

$$
\begin{align*}
E_1 + E_2 + \ldots + E_n &= \frac{\hat{E_1} + \hat{E_2} + \ldots + \hat{E_n} - CI}{1+C} 
\end{align*}
$$

So then for each $E_i$, we can calculate the true $E_i$ as:

$$
E_i = \frac{\hat{E_i}}{\sum_{j=1}^n \hat{E_j}} \times \frac{\sum_{j=1}^n \hat{E_j} - CI}{1+C}
$$




