# GTF_manager
My personal gtf manager module.

Iterate over gtf file and extract relevant informations
## DEPENDENCIES

- pybedtools
- tqdm

## BASIC USAGE

GTF_managers offers a easy-use parser to handle gtf files.

With simple commands you can creates a object called gtf and parses a file inside it.
```
    from GTF_manager import GTF_manager
    
    gtf = GTF_manager('path/to/my_gtf_file.gtf')
```

GTF_manager separates the file lines in two main object types: 'genes' and 'transcripts'.
The 'genes' contains the description about total locus and all keep a list of all transcripts inside that locus.
To print all gene locus in bed6 format in your screen you need the function `gene_list_to_bed6`:
```
    for genebed6 in gtf.gene_list_to_bed6(): #
       print genebed6.gene_id
```
Saves a bed6 files is easy:
```
    gtf.gene_list_to_bed6(file_name='path/to/save/my_experimental_bed6.gtf', save_in_file=True) # save bed6 in file
```
There the possibility to get a specific gene using the function `get_gene()`:

```
    my_gene = gtf.get_gene("ENSG00000180198.11") # Pass the gene_id with argument
```
Now `my_gene` objects is a instance of Gene_content class and you can do some operations.
Here are a small list of possible operations using this class.

``` 
    print type(my_gene)
    
    <type 'int'> #check the type
    
    print my_gene.locus_coords(): show the coords of locus
    print my_gene.locus_length() # length of locus
    print my_gene.get_bed6() # generate a bed6 format string
    print my_gene.number_of_transcripts # print the number of transcripts inside locus
 ```
 
 Getting the name of all transcripts inside locus:
 ```
    get_transcripts_ids()
    
```
You can iterate over all transcripts inside a specific locus and print yours `transcript_name` using :
```
    for transcripts in my_gene.get_transcripts():
        print transcripts.transcript_name
```

If you are interested in check the fields inside the gtf object, you can print this keys in the screen using the  'gtf' object.
```
    gtf.print_attrs_fields()  # list all attr possible inside Transcripts objects.
```
After gets all possible attrs you can call this attr inside `transcript` objects

```
    for transcripts in my_gene.get_transcripts():
        print transcripts.attrs("transcript_status")
```
    

