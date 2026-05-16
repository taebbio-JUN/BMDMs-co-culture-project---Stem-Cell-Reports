# Data

Raw and processed sequencing data for this project are **not stored in this repository** due to file-size limits. They are deposited in the appropriate public archives.

## Raw FASTQ files

- **Repository:** ArrayExpress (ENA mirror)
- **ArrayExpress accession:** `E-MTAB-XXXXX`
- **ENA Study accession:** `PRJEB-XXXXX`
- **URL:** https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-XXXXX

## Processed count matrix

- **File:** `gene_count_matrix.csv`
- **Distribution:** included as processed data with the ArrayExpress submission above
- **Additional archive:** Zenodo (DOI: `10.5281/zenodo.XXXXXXX`)

## Sample metadata

| Sample ID | Group | Replicate | Index 1 | Index 2 |
| --- | --- | :---: | --- | --- |
| Con1R | Organoid | 1 | ATCTATGG | TCCGTCAA |
| Con2R | Organoid | 2 | TCGAATAC | GCCAGGTG |
| Con3R | Organoid | 3 | ACATGGAG | AGACTCCT |
| 1K1R | Organoid_BMDM_1K | 1 | CCAGCTAT | CCGCTGTC |
| 1K2R | Organoid_BMDM_1K | 2 | ATGGACCT | GATAACAC |
| 1K3R | Organoid_BMDM_1K | 3 | CCACGGTC | CTAGACGC |
| 5K1R | Organoid_BMDM_5K | 1 | AGACGTTA | CGCTGCAT |
| 5K2R | Organoid_BMDM_5K | 2 | GACAATAA | CGCTCATC |
| 5K3R | Organoid_BMDM_5K | 3 | AGGCCTGG | AACACGTC |
| TNFa2 | Organoid_TNFa | 1 | TTAGCAGG | GTTACACA |
| TNFa3 | Organoid_TNFa | 2 | TACGCAAG | AACGGAAG |
| TNFa4 | Organoid_TNFa | 3 | CGTGGAGG | ATTACGTG |
| TNFa1K1 | Organoid_TNFa_BMDM_1K | 1 | TTCGTGCA | TGATATGA |
| TNFa1K2 | Organoid_TNFa_BMDM_1K | 2 | TGTGATCG | TTACTTAT |
| TNFa1K3 | Organoid_TNFa_BMDM_1K | 3 | TGACGTTG | AGTAGGCC |
| TNFa5K1 | Organoid_TNFa_BMDM_5K | 1 | AATGGTTA | TAGATGTT |
| TNFa5K2 | Organoid_TNFa_BMDM_5K | 2 | GAGAGACA | GAATGGTT |
| TNFa5K3 | Organoid_TNFa_BMDM_5K | 3 | TCGAGACA | GTTATCGA |
| BMDM1 | Macrophage | 1 | CAGGCAAC | GCAGAGTT |
| BMDM2 | Macrophage | 2 | TGTACTAA | CGTATGCA |
| BMDM3 | Macrophage | 3 | CGACGACG | TGTTCCTA |
| TNFaBMDM1 | Macrophage_TNFa | 1 | AGAATTAC | CTATAGCG |
| TNFaBMDM2 | Macrophage_TNFa | 2 | TCGGACTG | CGCAGTCG |
| TNFaBMDM3 | Macrophage_TNFa | 3 | CATCGGCC | TTGGATAC |

## Sequencing parameters

| Parameter | Value |
| --- | --- |
| Sequencing provider | Macrogen Inc. (Seoul, Republic of Korea) |
| Library kit | TruSeq Stranded Total RNA with Ribo-Zero H/M/R_Gold |
| Application | Transcriptome Resequencing |
| Mode | Paired-end, 101 bp × 2 |
| Run type | Throughput Run |
| Target output per sample | ~6 Gb |
| Reference genome | GRCm38 (mouse) |

## Reproducing the analysis

1. Download `gene_count_matrix.csv` from one of the archives listed above.
2. Place it in the project root directory.
3. Run the pipeline as described in the top-level `README.md`.
