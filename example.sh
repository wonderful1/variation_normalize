gzip -d test_data/chrY.fa.gz
python variation_normalize.py -ref test_data/chrY.fa -f test_data/chrY.vcf.gz -c 1,2,4,5  -o ./test.out
