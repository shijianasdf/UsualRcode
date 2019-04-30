suppressPackageStartupMessages({
  # Bioconductor
  library(BiocParallel)
  library(SingleCellExperiment)
  library(clusterExperiment)
  library(scone)
  library(zinbwave)
  library(slingshot)
  # CRAN
  library(gam)
  library(RColorBrewer)
})
#> Warning: replacing previous import 'SingleCellExperiment::weights' by
#> 'stats::weights' when loading 'slingshot'
#> Warning in rgl.init(initValue, onlyNULL): RGL: unable to open X11 display
#> Warning: 'rgl_init' failed, running with rgl.useNULL = TRUE
set.seed(20)
register(SerialParam())
library(fletcher2017data)

data(fletcher)
fletcher
#> class: SingleCellExperiment 
#> dim: 28284 849 
#> metadata(0):
#> assays(1): counts
#> rownames(28284): Xkr4 LOC102640625 ... Ggcx.1 eGFP
#> rowData names(0):
#> colnames(849): OEP01_N706_S501 OEP01_N701_S501 ... OEL23_N704_S503
#>   OEL23_N703_S502
#> colData names(19): Experiment Batch ... CreER ERCC_reads
#> reducedDimNames(0):
#> spikeNames(0):
colData(fletcher)
#> DataFrame with 849 rows and 19 columns
#>                          Experiment    Batch publishedClusters    NREADS
#>                            <factor> <factor>         <numeric> <numeric>
#> OEP01_N706_S501     K5ERRY_UI_96HPT      Y01                 1   3313260
#> OEP01_N701_S501     K5ERRY_UI_96HPT      Y01                 1   2902430
#> OEP01_N707_S507     K5ERRY_UI_96HPT      Y01                 1   2307940
#> OEP01_N705_S501     K5ERRY_UI_96HPT      Y01                 1   3337400
#> OEP01_N704_S507     K5ERRY_UI_96HPT      Y01                -2    117892
#> ...                             ...      ...               ...       ...
#> OEL23_N704_S510 K5ERP63CKO_UI_14DPT      P14                -2   2407440
#> OEL23_N705_S502 K5ERP63CKO_UI_14DPT      P14                -2   2308940
#> OEL23_N706_S502 K5ERP63CKO_UI_14DPT      P14                12   2215640
#> OEL23_N704_S503 K5ERP63CKO_UI_14DPT      P14                12   2673790
#> OEL23_N703_S502 K5ERP63CKO_UI_14DPT      P14                 7   2450320
#>                  NALIGNED    RALIGN TOTAL_DUP    PRIMER
#>                 <numeric> <numeric> <numeric> <numeric>
#> OEP01_N706_S501   3167600   95.6035   47.9943 0.0154566
#> OEP01_N701_S501   2757790   95.0167    45.015 0.0182066
#> OEP01_N707_S507   2178350   94.3852   43.7832 0.0219196
#> OEP01_N705_S501   3183720   95.3952   43.2688 0.0183041
#> OEP01_N704_S507     98628   83.6596   18.0576 0.0623744
#> ...                   ...       ...       ...       ...
#> OEL23_N704_S510   2305060   95.7472   47.1489 0.0159111
#> OEL23_N705_S502   2203300   95.4244   62.5638 0.0195812
#> OEL23_N706_S502   2108490   95.1637   50.6643 0.0182207
#> OEL23_N704_S503   2568300   96.0546   60.5481 0.0155611
#> OEL23_N703_S502   2363500   96.4567   48.4164 0.0122704
#>                 PCT_RIBOSOMAL_BASES PCT_CODING_BASES PCT_UTR_BASES
#>                           <numeric>        <numeric>     <numeric>
#> OEP01_N706_S501               2e-06          0.20013      0.230654
#> OEP01_N701_S501                   0         0.182461       0.20181
#> OEP01_N707_S507                   0         0.152627      0.207897
#> OEP01_N705_S501               2e-06         0.169514      0.207342
#> OEP01_N704_S507             1.4e-05         0.110724      0.199174
#> ...                             ...              ...           ...
#> OEL23_N704_S510                   0         0.287346      0.314104
#> OEL23_N705_S502                   0         0.337264      0.297077
#> OEL23_N706_S502               7e-06         0.244333      0.262663
#> OEL23_N704_S503                   0         0.343203      0.338217
#> OEL23_N703_S502               8e-06         0.259367      0.238239
#>                 PCT_INTRONIC_BASES PCT_INTERGENIC_BASES PCT_MRNA_BASES
#>                          <numeric>            <numeric>      <numeric>
#> OEP01_N706_S501           0.404205             0.165009       0.430784
#> OEP01_N701_S501           0.465702             0.150027       0.384271
#> OEP01_N707_S507           0.511416              0.12806       0.360524
#> OEP01_N705_S501           0.457556             0.165586       0.376856
#> OEP01_N704_S507           0.489514             0.200573       0.309898
#> ...                            ...                  ...            ...
#> OEL23_N704_S510           0.250658             0.147892        0.60145
#> OEL23_N705_S502           0.230214             0.135445       0.634341
#> OEL23_N706_S502           0.355899             0.137097       0.506997
#> OEL23_N704_S503           0.174696             0.143885        0.68142
#> OEL23_N703_S502           0.376091             0.126294       0.497606
#>                 MEDIAN_CV_COVERAGE MEDIAN_5PRIME_BIAS MEDIAN_3PRIME_BIAS
#>                          <numeric>          <numeric>          <numeric>
#> OEP01_N706_S501           0.843857           0.061028           0.521079
#> OEP01_N701_S501            0.91437            0.03335           0.373993
#> OEP01_N707_S507           0.955405           0.014606            0.49123
#> OEP01_N705_S501            0.81663           0.101798           0.525238
#> OEP01_N704_S507            1.19978                  0           0.706512
#> ...                            ...                ...                ...
#> OEL23_N704_S510           0.698455           0.198224           0.419745
#> OEL23_N705_S502           0.830816           0.105091           0.398755
#> OEL23_N706_S502           0.805627           0.103363           0.431862
#> OEL23_N704_S503           0.745201           0.118615            0.38422
#> OEL23_N703_S502           0.711685           0.196725           0.377926
#>                     CreER ERCC_reads
#>                 <numeric>  <numeric>
#> OEP01_N706_S501         1      10516
#> OEP01_N701_S501      3022       9331
#> OEP01_N707_S507      2329       7386
#> OEP01_N705_S501       717       6387
#> OEP01_N704_S507        60        992
#> ...                   ...        ...
#> OEL23_N704_S510       659          0
#> OEL23_N705_S502      1552          0
#> OEL23_N706_S502         0          0
#> OEL23_N704_S503         0          0
#> OEL23_N703_S502      2222          0