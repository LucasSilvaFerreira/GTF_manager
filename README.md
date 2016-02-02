# GTF_manager
My personal gtf manager module.

Iterate over gtf file and extract relevant informations
## DEPENDENCIES
- pybedtools
- tqdm
##BASIC USAGE

GTF_managers offers a easy-use parser to handle gtf files.
         Examples:
             gtf = GTF_manager('path/to/my_gtf_file.gtf')
             for genebed6 in gtf.gene_list_to_bed6(): #
                print genebed6.gene_id
             gtf.gene_list_to_bed6(file_name='path/to/save/my_experimental_bed6.gtf', save_in_file=True) # save bed6 in file
             gtf.gene_list_to_bed6()
             gtf.print_attrs_fields()  # list all attr possible inside Transcripts objects.
             my_gene = gtf.get_gene("ENSG00000180198.11") # use gene_id to get a Gene_content object.
             print my_gene.locus_length() # length of locus
             print my_gene.transcripts_list() # List gene transcripts
             print my_gene.get_bed6() # generate a bed6 format string
             for transcripts in my_gene.get_transcripts():
                 print transcripts.attrs("transcript_status")
