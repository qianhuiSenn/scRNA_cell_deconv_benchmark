# Dataset
The four pairs of Real datasets are available through: https://drive.google.com/drive/folders/1famn_tRVT_Es0GgFng508QggXcqETPCJ?usp=sharing
# Simulation Datasets
Please find the corresponding Rscripts for the generation of the simulation data.
# Installation of the Docker Files
Pull the software version used in the manuscript in an R session.
```bash
docker pull qhhuang/benchmark_celltype_r_packages
mkdir /<Results data folder>/
cd /<Results data folder>/
PATHDATA=`pwd`
```
## usage

The pipeline consists of 3 steps (for downloading the data) and 4 steps for aligning and calling SNVs:

```bash
# Download
docker run --rm opoirion/ssrge download_soft_file -h
docker run --rm opoirion/ssrge download_sra -h
docker run --rm opoirion/ssrge extract_sra -h
# align and SNV calling
docker run --rm opoirion/ssrge star_index -h
docker run --rm opoirion/ssrge process_star -h
docker run --rm opoirion/ssrge feature_counts -h
docker run --rm opoirion/ssrge process_snv -h

```

## example

Let's download and process 2 samples from GSE79457 in a project name test_n2

```bash
# download of the soft file containing the metadata for GSE79457
docker run --rm -v $PATHDATA:/data/results/:Z opoirion/ssrge download_soft_file -project_name test_n2 -soft_id GSE79457
# download sra files
docker run --rm -v $PATHDATA:/data/results/:Z opoirion/ssrge download_sra -project_name test_n2 -max_nb_samples 2
# exctract sra files
docker run --rm -v $PATHDATA:/data/results/:Z opoirion/ssrge extract_sra -project_name test_n2
# rm sra files (optionnal)
docker run --rm -v $PATHDATA:/data/results/:Z opoirion/ssrge rm_sra -project_name test_n2
## all these data can also be obtained using other alternative workflows
# here you need to precise which read length to use for creating a STAR index and which ref organism (MOUSE/HUMAN)
docker run --rm -v $PATHDATA:/data/results/:Z opoirion/ssrge star_index -project_name test_n2 -read_length 100 -cell_type HUMAN
# STAR alignment
docker run --rm -v $PATHDATA:/data/results/:Z opoirion/ssrge process_star -project_name test_n2 -read_length 100 -cell_type HUMAN
# sample-> gene count matrix
docker run --rm -v $PATHDATA:/data/results/:Z opoirion/ssrge feature_counts -project_name test_n2
#SNV inference
docker run --rm -v $PATHDATA:/data/results/:Z opoirion/ssrge process_snv -project_name test_n2 -cell_type HUMAN
```
