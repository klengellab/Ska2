**SKA2 regulated hyperactive secretory autophagy drives neuroinflammation induced neurodegeneration**


**Abstract**

High levels of proinflammatory cytokines induce neurotoxicity and catalyze inflammation-driven neurodegeneration, 
but the specific release mechanisms from microglia remain elusive. We demonstrate that secretory autophagy (SA), 
a non-lytic modality of autophagy for secretion of vesicular cargo, regulates neuroinflammation-mediated neurodegeneration via SKA2 and FKBP5 signaling. 
SKA2 inhibits SA-dependent IL-1β release by counteracting FKBP5 function. Hippocampal Ska2 knockdown in male mice hyperactivates SA resulting in neuroinflammation, 
subsequent neurodegeneration and complete hippocampal atrophy within six weeks. 
The hyperactivation of SA increases IL-1β release, contributing to an inflammatory feed-forward vicious cycle including NLRP3-inflammasome activation and Gasdermin D-mediated neurotoxicity, 
which ultimately drives neurodegeneration. Results from protein expression and co-immunoprecipitation analyses of male and female postmortem human brains demonstrate that SA is hyperactivated in Alzheimer’s disease. Overall, our findings suggest that SKA2-regulated, hyperactive SA facilitates neuroinflammation and is linked to Alzheimer’s disease, providing new mechanistic insight into the biology of neuroinflammation.

**Data and analysis**

* Raw data: RNAseq fastq files can be accessed through GEO accession GSE181203
* Processed data: The data that is used for the scripts is uploaded to the github page, the counts table is stored in the counts.csv, metadata is stored in the metadata.csv. Additionally, we also                         provide the output from MuSiC were the cell type proportions are stored in the all_cell_type_prop.csv.
* Analysis: The main analysis normalization, differential expression analysis and enrichment analysis is performed in the RNAseq_analayis.R
            and the deconvolution of the RNAseq data with MuSiC is performed in the Deconvolution_MuSiC.R. Additionally, for the reviewer we also
            explored other Deseq2 models and the code is available for those models in RNAseq_analysis_multiple_DEseq2_models.R. The functions that are used throughout the analysis can be found in the               functions.R script
  
