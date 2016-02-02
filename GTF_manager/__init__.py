__author__ = 'lucas'
# -*- coding: utf-8 -*-
import sys
import re
from tqdm import tqdm
from pybedtools import BedTool, tempfiles
import cPickle

class Transcript():
    def __init__(self, transcript_name, exons_array, attr_list):
        self.transcript_name = transcript_name
        self.exons = exons_array
        self.gene_info = self.exons[0][1]
        self.strand = self.exons[0][6]
        self.__attr_list = attr_list
        self.exon_count = len(exons_array)
        self.__attrs = self.__parse_attrs()  # aqui eu deixarei todos os atributos dos genes em um dict



    def get_attrs(self):
        return self.__attrs

    def set_new_attr(self, attr_key, attr_value):
        dict_attr = self.__attrs
        dict_attr[attr_key] = attr_value

    def locus_size_total(self):
        chr, s, e = self.locus_coords()
        return e-s

    def transcript_exons_length_sum(self):
        return self.__get_transcrip_size()

    def exon_size(self, one_exon_in_array):
        return int(one_exon_in_array[4]) - int(one_exon_in_array[3])

    def __get_transcrip_size(self):
        return sum([self.exon_size(one_exon_in_array=exon) for exon in self.exons])

    def __parse_attrs(self):

        # selecting just the unique fields inside all transcripts exons.
        fields_unique={}
        for field_u_exon in self.exons:
            field_attrs = field_u_exon[-1].strip(';').split(';')
            for field_x in field_attrs:
                field_key = field_x.strip(' ').split(' ')[0]
                field_value = field_x.strip(' ').split(' ')[1]
                if field_key in fields_unique:
                    fields_unique[field_key].append(field_value)
                else:
                    fields_unique[field_key] = []
                    fields_unique[field_key].append(field_value)
        #poping the keys with size different of exons size array. THis will keep just no variant attrs

        keys_to_pop = [possible_key for possible_key, value in fields_unique.iteritems() if len(list(set(value)))!=1]

        #exclude keys_to_pop from attrs
        fields_to_parse= {}
        for field in self.exons[-1][-1].strip(';').split(';'):
            key_parse = field.strip(' ').split(' ')[0]
            value_parse = field.strip(' ').split(' ')[1]
            if not key_parse in keys_to_pop:
                fields_to_parse[key_parse] = value_parse
            else: # this happen when the key change between exons
                fields_to_parse[key_parse] = "Non_transcript_level_information"


        hash_to_load_return = {}
        for attr in self.__attr_list:
            if attr in fields_to_parse:
                hash_to_load_return[attr] = fields_to_parse[attr]
            else:
                hash_to_load_return[attr] = None

        return hash_to_load_return

    def locus_coords(self):
        chr = self.exons[0][0]
        strand = self.exons[0][6]

        if strand == "-":

            return (chr,
                    min([int(end[3]) for end in self.exons]),
                    max([int(start[4]) for start in self.exons]))
        elif strand == "+":
            return (chr,
                    min([int(start[3]) for start in self.exons]),
                    max([int(end[4]) for end in self.exons]))

    def get_bed6(self):
        chr, start, end = map(str, self.locus_coords())
        return '\t'.join([chr, start, end, self.transcript_name, self.gene_info, self.strand])

class Gene_content():
    def __init__(self, gene_id,
                 strand,
                 transcripts_id_hash,
                 info_gene_str,
                 attr_transcript_field
                 ):
        self.gene_id = gene_id
        self.gene_info = info_gene_str
        self.strand = strand  # checar se strand existe na posicao correta else: return error
        self.transcripts_ids = transcripts_id_hash  # sorted_dict
        self.__child_transcript_possible_fields = attr_transcript_field  # hash com os valores que os transcritos vao possuir
        self.__transcripts_list = self.__parse_transcripts()
        self.number_of_transcripts = len(self.__transcripts_list)

    def get_bed6(self):
        chr, start, end = map(str, self.locus_coords())
        return '\t'.join([chr, start, end, self.gene_id, self.gene_info, self.strand])

    def locus_coords(self):
        # print type(self.__transcripts_list), self.__transcripts_list
        chr = self.__transcripts_list[0].exons[0][0]
        strand = self.__transcripts_list[0].exons[0][6]
        if strand == "-":

            return (chr,
                    min([int(end[3]) for end_parse in self.__transcripts_list for end in end_parse.exons]),
                    max([int(start[4]) for start_parse in self.__transcripts_list for start in start_parse.exons]))
        elif strand == "+":
            return (chr,
                    min([int(start[3]) for start_parse in self.__transcripts_list for start in start_parse.exons]),
                    max([int(end[4]) for end_parse in self.__transcripts_list for end in end_parse.exons]))

    def locus_length(self):
        chr, start_l, end_l = self.locus_coords()
        return end_l - start_l

    def transcripts_list(self):
        return self.__parse_transcripts()

    def get_transcripts_ids(self):
        out_ids_trans = []

        for trans_line in self.transcripts_ids.iterkeys():
            # print trans_line
            out_ids_trans.append(re.search('transcript_id "(\S*)"', trans_line[-1]).group(1))
        return out_ids_trans

    def __parse_transcripts(self):
        hash_transcripts_to_parse = {}
        return [Transcript(transcript_name=ts_key,
                           exons_array=ts_values_array,
                           attr_list=self.__child_transcript_possible_fields) for ts_key, ts_values_array in
                self.transcripts_ids.iteritems()]

class GTF_manager():
    def __init__(self, gtf_file):
        """GTF_managers offers a easy-use parser to handle gtf files.
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

        Attributes:
            self.file_name (str): Return the GTF file PATH


        Args:
            gtf_file (str): system PATH with a gtf file.
        """
        print 'O PARSER DE ATTRS deve mostrar apenas os campos que se repetem. Ou os attrs ficariam com informacoes que modificam entre os de exons'
        self.file_name = gtf_file
        self.gtf_file_open = [line_gtf.split('\t') for line_gtf in open(gtf_file, 'r').read().split('\n') if
                              len(line_gtf) > 0]
        self.__attr_list = {} # pode dar erro depois de extrair genes utilizando select
        self.genes_hash = self.__parse_lines()

    def __parse_lines(self):
        hash_gtf_file = {}
        print 'Parsing genes...'

        for line in tqdm(self.gtf_file_open):
            if line:
                # Populating attrs list
                for attr_line in line[-1].strip(';').split(
                        ';'):  # retirar ultima ocorrencia de ponto-virgula e splitar por ponto-virgula
                    self.__attr_list[attr_line.strip().split(" ")[0]] = ''
                # Parsing lines
                try:
                    gtf_gene_id = re.search('gene_id "(\S*)"', line[-1]).group(1)
                    if gtf_gene_id in hash_gtf_file:
                        hash_gtf_file[gtf_gene_id].append(line)
                    else:
                        hash_gtf_file[gtf_gene_id] = []
                        hash_gtf_file[gtf_gene_id].append(line)
                except:
                    raise IOError("\n The Line: \n{line}\n don't have a regular gene_id field!".format(line=line))

        gene_list_to_return = {}
        print '\n'
        print 'Converting gene to proper format...'
        for gene_unpack in tqdm(hash_gtf_file):
            key_gene = gene_unpack
            gene_value = hash_gtf_file[key_gene]
            gene_id = key_gene
            gene_info = gene_value[0][1]
            strand = gene_value[0][6]
            transcript_hash = {}
            for transcript in gene_value:  # Loading  each exon in each transcript
                id_line = re.search('transcript_id "(\S*)"', transcript[-1]).group(1)

                if id_line in transcript_hash:
                    transcript_hash[id_line].append(transcript)
                else:
                    transcript_hash[id_line] = []
                    transcript_hash[id_line].append(transcript)
            # Gene content class creation
            gene_object = Gene_content(gene_id=gene_id,
                                       strand=strand,
                                       transcripts_id_hash=transcript_hash,
                                       info_gene_str=gene_info,
                                       attr_transcript_field=self.__attr_list)

            if gene_id in gene_list_to_return:
                raise IOError('\n Duplicated gene_id\n {gene_id}\n'.format(gene_id=gene_id))
            else:
                gene_list_to_return[gene_id] = gene_object

        return gene_list_to_return

        print '\n Parsing finished...'

    def number_of_genes_loci(self):
        return len(self.genes_hash)

    def gene_list(self):
        """

        Returns:
            Gene_content
        """
        return [gene_content for gene_content in self.genes_hash.itervalues()]

    def gene_list_to_bed6(self, file_name=None, save_in_file=False):
        """

        Returns:

            BedTool: BedTool object with gene_ids inside.
            or
            Save a ".bed" in specified path in save_in_file parameter.
        """
        bed6_genes = [gene_to_b6.get_bed6() for gene_to_b6 in self.gene_list()]
        # aqui deve entrar algum distema de filtros
        bed6_genes = BedTool('\n'.join(bed6_genes), from_string=True).sort()

        if not save_in_file:
            return bed6_genes
        else:
            if file_name:
                return bed6_genes.saveas(fn=file_name, trackline="track name='Genes {}' color=128,0,0".format(
                        file_name.split('/')[-1]))
            else:
                raise IOError('\nthe file_name method function needs a name or complete file path with name.\n')

    def gene_list_select(self, arr_gene_id_to_select):
        ids_not_fount = []
        hash_to_modify = {}
        for remove_id in arr_gene_id_to_select:
            try:
                hash_to_modify[remove_id] = self.genes_hash[remove_id]
            except:
                ids_not_fount.append(remove_id)
                print remove_id
        if len(ids_not_fount) > 0:
            print len(ids_not_fount), 'Ids não encontrados...'
        copy_to_return = cPickle.loads(cPickle.dumps(self, -1))
        copy_to_return.genes_hash = hash_to_modify
        return copy_to_return

    def transcripts_list(self):
        transcript_list_return = []
        for gene_in_list in self.gene_list():
            for transcript_in_list in gene_in_list.transcripts_list():
                transcript_list_return.append(transcript_in_list)
        return transcript_list_return

    def transcripts_list_to_bed6(self, file_name=None, save_in_file=False):
        """

        Returns:

            BedTool: BedTool object with Transcripts inside.
            or
            Save a ".bed" in specified path in save_in_file parameter.
        """

        bed6_trans = [trans_to_b6.get_bed6() for trans_to_b6 in self.transcripts_list()]
        # aqui deve entrar algum distema de filtros
        bed6_trans = BedTool('\n'.join(bed6_trans), from_string=True).sort()

        if not save_in_file:
            return bed6_trans
        else:
            if file_name:
                return bed6_trans.saveas(fn=file_name, trackline="track name='Transcripts {}' color=128,0,0".format(
                        file_name.split('/')[-1]))
            else:
                raise IOError('\nthe file_name method function needs a name or complete file path with name.\n')

    def get_gene(self, gene_id):
        return self.genes_hash[gene_id]


    def print_attrs_fields(self):

        for attr_field_print in self.__attr_list.iterkeys():
            print attr_field_print
