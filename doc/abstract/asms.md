
## Title – Limit 20 words

Triqler for Data Independent Aquisition Data

## Introduction – Limit 120 words


120 words

Shotgun proteomics is frequently used for determining differentially abundant proteins between experimental conditions.There are two common acquisition methods in shotgun proteomics, data-dependent acquisition (DDA), and data-independent acquisition (DIA). In DDA, the mass spectrometer selects peptide ions to be fragmented and analyzed based on intensity, while in DIA, the instrument is scheduled to fragment larger mass windows in an unbiased manner.

We have previously designed a hierarchical Bayesian model, Triqler, able to control for errors from both the identification and quantification process in LFQ experiments. Triqler integrates the error probabilities from identification and quantification, to obtain better accuracy in calling differentially abundant proteins.  However, Triqler was developed with DDA-data in mind. Here we demonstrated that Triqler also is compatible with DIA-data.

## Methods – Limit 120 words

126 words

Triqler assumes that peptide abundances follow a probability distribution that is initiated according to a prior distribution, but updated based on each sample's registered peptide abundance values. This results in peptide abundance distributions that are integrated into protein abundance and subsequently fold change distributions. In the processing Triqler weights in information on differences in peptide abundances within sample groups and search engine identification error probabilities. Data that indicate uncertainty results in wider abundance distributions, while certainty results in tighter abundance distributions. Triqler integrates the resulting fold change distributions into posterior probabilities and q-values for each protein having an actual fold change larger than a preset threshold value.

## Preliminary data – Limit 300 words

Using the data from LFQBench we showed that DIA data conforms with Triqler assumptions that the structure of measurements is dominated by multiplicative errors and that missing peptide abundance follows a censored normal distribution and that we hence can use the method also for DIA data.

We further benchmarked Triqler against a set of state-of-the-art methods for protein summarization methods, MSstats, MSqRob2, and Top3. Comparisons were made both on sets with known concentrations i.e. LFQBench, and some cancer proteomics-related datasets. To reduce the effect of the different protein inference strategies for the tested protein summarization tools, a modified FASTA file, without shared peptides, was used for database search. The filter removed protein sequences with shared peptides so that the final database contained no two proteins sharing tryptic peptides longer than 7 amino acids.

For all the compared sets and processing options, Triqler was found to outperform the other protein summarization tools. I was also found to be statistically well-calibrated. When using the LFQBench dataset with known relative protein abundances, Triqler's reported quantification errors that corresponded better to the actual error rates than the compared methods.

## Novel aspect – Limit 20 words

Triqler is compatible with DIA data. 




To ensure a fair comparison we used two pipelines for processing the data and reduced the effect of protein inference by removing proteins with shared peptides.  

We compared Triqler to other state-of-the-art summarization methods and found that Triqler was able to detect a higher number of differentially abundant proteins and made fewer errors, see Figure 1. Further, we argue that Triqler is both easier to use and induces less bias as it does not require any filtering or imputation before analysis.

As proteomics is complex, computational methods that are user-friendly and reproducible are important. We demonstrate that Triqler improves the field of quantitative proteomics in both regards.
 
