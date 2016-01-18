__author__ = 'lucas'
# -*- coding: utf-8 -*-
import sys
from GTF_manager import GTF_manager
import cPickle







def main():
    arquivo = GTF_manager("/work/users/dinar/genome_human/Ensembl_GRCh37/genes.gtf")
    print 'number', arquivo.number_of_genes_loci()
    my_copy = arquivo.gene_list_select(['ENSG00000180198'])
    print 'number', my_copy.number_of_genes_loci()
    for x in my_copy.gene_list():
        print x.gene_id

if __name__ == '__main__':
    sys.exit(main())
