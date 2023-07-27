# ARC_Technical_Exercise

## Manual First Commands

Install bbduk from sourceforge: https://sourceforge.net/projects/bbmap/

add to Path

```
export PATH="/workspace/bbmap_install/bbmap":$PATH
```

### Qiime2 & mOTUs

Install Qiime2 via anaconda as per https://docs.qiime2.org/2023.5/install/native/#install-qiime-2-within-a-conda-environment

```
wget https://data.qiime2.org/distro/core/qiime2-2023.5-py38-linux-conda.yml
mamba env create -n qiime2-2023.5 --file qiime2-2023.5-py38-linux-conda.yml
rm qiime2-2023.5-py38-linux-conda.yml # cleanup
```

Install q2-mOTUs as per https://github.com/motu-tool/q2-mOTUs#installation

```
git clone https://github.com/motu-tool/q2-mOTUs
cd q2-mOTUs
make install
```



### Deliverables

A csv or txt report including the following fields:

* Number of molecules sequenced {parse from fastqc}
* Sequencing configuration (read 1, index 1, index 2, and read 2 lengths) {indexes from direct fastq, lengths from fastqc}
* Library type (Nextera or TruSeq) {parse from fastqc}
* Mean and standard deviation of library insert sizes after adapter removal {bbmerge after bbduk}
* RefSeq or Genbank accession numbers of the sequenced organism. {entrez-direct approach}
  * Identify taxon ids via assembly-based approach & mmseq2 taxonomy
  * Compare to read based approach quimme2+motus
  * `esearch -db taxonomy -query "txid320324[Organism:exp]" | elink -target assem
bly | efetch -format docsum  | xtract -pattern DocumentSummary -element RefSeq`

---

* A plot showing the instrument reported phred scores and alignment inferred quality scores across all positions of read 1 and read 2. Read 1 and Read 2 can be combined on the same plot or separated as 2 plots. The x-axis should represent positions across the read in units of base pairs. The y-axis should represent quality scores, calculated as -10*log10(error rate). Two lines on each plot should represent mean phred qualities reported by the instrument and quality scores inferred from alignment error rates respectively.
  * Downloading genome fna for build:
    * `wget "$(esearch -db assembly -query "GCF_000716445.1" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F"/" '{print $0"/"$NF"_genomic.fna.gz"}')"`
  * Align with bowtie2 to top representative genomes from NCBI and plot
    * build: `bowtie2-build /workspace/entrez_direct/GCF_000716445.1_ASM71644v1_genomic.fna Streptomyces_wedmorensis`