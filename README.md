### variation_normalize
 a variant normalization procedure


### usage
View instructions by running `python3 variation_normalize.py -h`

```
python3 variation_normalize.py [-h] -ref FASTA -f info -c column number -o outpre

optional arguments:
  -h, --help            show this help message and exit
  -ref FASTA, --reference FASTA
                        The genome sequencing
  -f info, --input_file info
                        the path of input file, should be Tab-separated file.
  -c column number, --cpra column number
                        the column number of chr-pos-ref-alt(CPRA), for example: 1,2,4,5
  -o outpre, --outpre outpre
                        filename of output
```

### example
```BASH
python variation_normalize.py -ref test_data/chrY.fa -f test_data/chrY.vcf.gz -c 1,2,4,5  -o ./test.out

```
