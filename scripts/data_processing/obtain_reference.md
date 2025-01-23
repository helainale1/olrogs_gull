# Download Reference from NCBI 
https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_964199755.1/

`scp path/in/local username@cluster:path/in/cluster`
```bash
scp /Users/hdl5108/desktop/ncbi_dataset.zip hdl5108@submit.hpc.psu.edu:/storage/group/dut374/default/helaina/data/L_michahellis_ref

```
## unzipping 
since the folder is zipped, I will unzip it using the `unzip` command
```bash
unzip /storage/group/dut374/default/helaina/data/L_michahellis_ref/ncbi_dataset.zip

```
## moving fna file into new folder

```bash
cp GCA_964199755.1_bLarMic1.1_genomic.fna /storage/group/dut374/default/helaina/data/L_michahellis_ref

```

