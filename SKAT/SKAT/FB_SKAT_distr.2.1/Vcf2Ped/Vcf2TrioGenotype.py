#!/usr/bin/env python

###############################################################################
###############################################################################
#  Last Updated Mar 27, 2013
#  Authors: Vlad Makarov
#  Language: Python
#  OS: UNIX/Linux, MAC OSX

## FREE FOR ACADEMIC AND OTHER NON-COMMERCIAL USERS

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.



###############################################################################
################################################################################

##   CONFIGURATION SECTION, MODIFY AS EXPLAINED BELOW IF NEEDED

###############################  KEYS (fields names) to MODIFY     #############

# AnnTools Keys used for the example attached
#KEYS=['name2', 'functionalClass', 'spliceInfo', 'spliceDist']

# SnpEff keys
KEYS=['SNPEFF_GENE_NAME', 'SNPEFF_FUNCTIONAL_CLASS', 'SNPEFF_EFFECT',  'SNPEFF_IMPACT']

# Annovar keys
#KEY=['subtype', 'transcript', 'classification', 'gene', 'genomechange', 'proteinchange']

# Keys are used by "parse_field(text, key, sep1, sep2)" method to separate fields:
# Ex name=NM_148901;name2=TNFRSF18;transcriptStrand=-;positionType=CDS;

SEP1=';'
SEP2='='

# Unknown genotypes ./.

UNKNOWN='0'

## File extensions used by this application
## It is not necessary to change them, unless there is a valid reason
## GENECOL='SNPEFF_GENE_NAME' #SnpEff
## GENECOL='gene' #Annovar

GENECOL='SNPEFF_GENE_NAME' #AnnTools
GENOTYPE_PED_EXT='ped'
VAR_EXT='var'
GENE_EXT='gene'
GENE_EFF_COL='SNPEFF_GENE_NAME' 



#############################   NO NEED TO MODIFY BELOW THIS POINT #############
###############################################################################
###############################################################################

import sys
import os.path
import os
import getopt
import gzip
from time import time



def print_help():
    h="""

    Purpose:  1. Creates GT (genotype) file based on the VCF and PED files provided
              for each trio listed in PED file
              2. Generates VARIANT file based on annotated VCF provided. By default, the snpEff/GATK
              annotation is expected, but can be configured to use other annotator. See CONFIGURATION SECTION
              for details
              3. Generated gene index file, i.e. start and end position of gene in the VCF file. Both positions are inclusive and start with one.


    Input:  1. Formatted and sorted variant call format (VCF) file. (http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41)
                Use VCFTools (http://vcftools.sourceforge.net/) for format verification if in doubt
                VCF file may be compressed.

            2. PED file. Must have the following columns IN THAT PARTICULAR ORDER:
                 FamilyID (0=unknown)
                 IndividualID (Must include all in the VCF file)
                 PaternalID (0=unknown)
                 MaternalID (0=unknown)
                 Sex (1=male; 2=female; 0=unknown)
                 Phenotype (1=control/unaffected, 2=case/affected, 0=unknown)

                PED format is described at http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped
            3. Base (no extension) of the output file

            4. True for PASS only, False for All variants

            5. True for minus 1 from Phenotype, as required by FB-SKAT

    Output: For each complete trio (child, father and mother ID present in both PED and VCF files),
    program outputs one row, where columns are:
    FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype, followed by Genotype1....GenotupeN
    First six columns are drawn from the PED file, the genotypes are drawn from matching samples
    at VCF file.

    Genotypes are coded as 0, 1 or 2 (no mutation, heterozygous, homozygous respectively)

    Note! Families with more than one child are treated as 2 or more trios.

    Dependencies:
    Python 2.6 or higher. No special libraries are needed

    Installation:
    No installation is needed. Copy the script to any directory at your computer and give it 755 permission

    For help:

    Bugs/suggestions may be sent to makarovv at gmail.com

    """
    print h


####################### General helper methods ################################
###############################################################################

""" Linear search, returns first index """
def find_first_index(lst, elem):
    if elem in lst:
        return lst.index(elem)
    return -1

def sublist(lst, index2include):
    slist=[]
    for i in index2include:
        slist.append(lst[i])
    return slist

def genotype_sum(lst):
    notZero=0
    for l in lst:
        if str(l).startswith('1/0') or str(l).startswith('0/1') or str(l).startswith('1/1') or str(l).startswith('1|0') or str(l).startswith('0|1') or str(l).startswith('1|1'):
            notZero=notZero+1
    return notZero

def minus_one(lst):
    newlst=[]
    for l in lst:
        newlst.append(str(int(l)-1));
    return newlst



def getFh(filename):
  fh = open(filename, "r")
  if filename.endswith('gz'):
    fh = gzip.open(filename, 'rb')
  return fh

""" Converts string to boolean """
def str2bool(v):
  return v.lower() in ["y", "yes", "true", "t", "1"]

""" Helper method to deduplicate the list while preserving the order"""
def dedup(mylist):
    outlist = []
    for element in mylist:
        if element not in outlist:
            outlist.append(element)
    return outlist

""" Reads specified column from the file"""
def read_one_str_col(filename, col=0, row=-1, sep='\t'):
    fh = open(filename, "r")
    values = []
    i=0
    for line in fh:
        if i>row:
            line=line.strip()
            if len(line)>0:
                fields = line.split(sep)
                values.append(fields[col].strip())
        i=i+1
    return values




""" Delete file """
def delete(filename):
    if os.path.exists(filename) and os.path.isfile(filename):
        os.unlink(filename)

""" read one column """
def col(filename, col, sep='\t'):
    fh = open(filename, "r")
    values = []
    for line in fh:
        line=line.strip()
        if len(line)>0:
            fields=line.split(sep)
            f=fields[col]
            values.append(f)
    return values



""" Replaces A,C,G,T to 1,2,3,4 whenever numeric values are needed"""
def toNumeric(allele):
    if allele=='A':
        return 1
    elif allele=='C':
        return 2
    elif allele=='G':
        return 3
    elif allele=='T':
        return 4
    else:
        return 0


###############################################################################
############################## Methods to genotype the trios ##################

""" transpose file """
def t(filein, fileout, sep='\t'):

    fh = open(filein)
    fh_out = open(fileout, "w")
    num_of_cols=len(fh.readline().strip().split(sep))
    print ("num_of_cols " + str(num_of_cols))
    for x in range(0,num_of_cols):
        colstr= sep.join(col(filename=filein, col=x, sep=sep))
        fh_out.write(colstr+'\n')
        print('Transposing ' + str(x+1) +' of ' +  str(num_of_cols))

""" make a list of samples in the VCF file"""
def make_list_of_samples(vcffile):
    fh=getFh(filename=vcffile)
    for line in fh:
        line = line.strip()
        if line.startswith('CHROM') or line.startswith('#CHROM'):
            fields = line.split('\t')
            samples=fields[9:len(fields)]
            print 'Vcf file contains ' + str(len(samples)) + ' samples'
            return samples



""" Loop through the PED file and return only families in order Child, Father, Mother"""
def make_list_of_trios(pedfile, vcffile):
    fh=getFh(filename=pedfile)
    include_samples_list=[]
    vcf_samples = make_list_of_samples(vcffile=vcffile)
    for line in fh:
        line = line.strip()
        if line.startswith('#')==False:
            #print line
            fields=line.split()
            IndividualID = str(fields[1]).strip()
            PaternalID = str(fields[2]).strip()
            MaternalID  = str(fields[3]).strip()
            #Sex = str(fields(4))
            if IndividualID != '0' and PaternalID != '0' and MaternalID != '0' and IndividualID in vcf_samples and PaternalID in vcf_samples and MaternalID in vcf_samples:
                include_samples_list.append(IndividualID)
                include_samples_list.append(PaternalID)
                include_samples_list.append(MaternalID)
            elif ( IndividualID != '0' and PaternalID != '0' and MaternalID != '0'):
                print IndividualID  +" "+  PaternalID +" "+ MaternalID;

    print 'PED file contains ' + str(len(include_samples_list)) + ' samples ('+ str(len(include_samples_list)/3) + ' trios)'
    return include_samples_list


""" Reformats the VCF file in order to:
        1. Include trios only
        2. Re-order samples in order Child, Father, Mother for each trio
"""
def reshufle_samples(filename, pedfile, passonly, comment_char='##'):

    include_samples_list=make_list_of_trios(vcffile=filename, pedfile=pedfile)
    #print include_samples_list
    tmp_vcf=filename+'.trios.tmp.vcf'
    fh_tmp_out = open(tmp_vcf, "w")

    print "Preparing tmp trio VCF file "
    print "Input file " + filename


    #include first 9 columns
    index2include=[0,1,2,3,4,5,6,7,8]

    fh=getFh(filename)
    for line in fh:
        line = line.strip()
        if line.startswith(comment_char):
            fh_tmp_out.write(line.strip()+'\n')


        elif line.startswith('CHROM') or line.startswith('#CHROM'):
            fields = line.split('\t')


            if len(fields)<10:
                print ('')
                print "######## ERROR PARSING FILE ##############"
                print "## VCF file format must have at least 10 fields: "
                print "## CHROM POS     ID        REF ALT    QUAL FILTER INFO   FORMAT"
                print "## followed by sample names"
                print "## See http://www.1000genomes.org/wiki/Analysis/vcf4.0 for documentation"
                print "## You may use our parse_info_table.py script to parse INFO in tabular file format "
                print "##########################################"
                break

            else:
                samples=fields[9:len(fields)]
                for s in include_samples_list:
                    i = int(find_first_index(samples, str(s).strip()))
                    if i > -1:
                        index2include.append(i+9)

                slist=sublist(fields, index2include)
                ln='\t'.join(slist).strip()
                fh_tmp_out.write(ln+'\n')




        else:
            fields = line.split('\t')
            ref=fields[3].strip()
            alt=fields[4].strip()
            #print "REF-ALT " + ref + '-' + alt
            if len(ref)==1 and len(alt)==1:
                slist=sublist(fields, index2include)
                genosum=genotype_sum(slist[9:len(slist)])
                filt=fields[6]

                if genosum > 0:
                    if passonly==True:
                        if filt=='PASS' or filt=='RefCoding':
                            ln='\t'.join(slist).strip()
                            fh_tmp_out.write(ln+'\n')
                    else:
                        ln='\t'.join(slist).strip()
                        fh_tmp_out.write(ln+'\n')


    fh.close()
    fh_tmp_out.close()


    print "Result is written to " + tmp_vcf


###>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

"""
    Generate standard ped file (2,1,0,9 for Homo, Het, Same as Ref, Unknown (./.))
    other specs are as http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped
"""


def vcf2pedSimple(vcf_file, pedfile, output_file_base, changePheno, header_line_start='#CHROM', comment_chars='##', sep='\t'):
    fh = open(vcf_file, "r")
    tmp_ped = output_file_base + '.tmp.onecol.'+GENOTYPE_PED_EXT
    fh_tmp_out = open(tmp_ped, "w")
    ped =  output_file_base + '.' + GENOTYPE_PED_EXT

    print "Generating genotyping ped file "
    print "Input file " + vcf_file


    FamilyIDlst = read_one_str_col(filename=pedfile, col=0)
    IndividualIDlst = read_one_str_col(filename=pedfile, col=1 )
    PaternalIDlst = read_one_str_col(filename=pedfile, col=2)
    MaternalIDlst = read_one_str_col(filename=pedfile, col=3)
    Sexlst = read_one_str_col(filename=pedfile, col=4)
    Phenotypelst = read_one_str_col(filename=pedfile, col=5)

    # Needed for FB-SKAT
    if changePheno==True:
        Phenotypelst=minus_one(Phenotypelst)


    required_fields=[FamilyIDlst, IndividualIDlst, PaternalIDlst, MaternalIDlst,Sexlst,Phenotypelst]

    for line in fh:
        ## skip comments
        if line.startswith(comment_chars):
            pass
        ## on the header line
        elif line.startswith(header_line_start) :
            #print line
            line = line.strip()
            fields = line.split(sep)
            samples=fields[9:len(fields)]
            #print samples
            inds=[]
            # make a list of samples and find an index of each in the sample info file
            # if sample is present in the VCF, but not in PED file, give warning
            for s in samples:
                ind = find_first_index(IndividualIDlst, s.strip() )
                if ind > -1:
                    inds.append(ind)
                else:
                    print ("Warning: Sample " + str(s) + " not found in the PED file " )



            # for each required field (FamilyID, PaternalID, MaternalID, IndividualID, Sex, Phenotype) , add values in order
            # of appearance in the VCF file
            for r in range(0, len(required_fields)):
                rflds=[]
                for i in inds:
                    rflds.append(required_fields[r][i])
                ## print them line by line to file
                fh_tmp_out.write(sep.join(rflds)+'\n')

        else:
            ## Genotype
            line = line.strip()
            fields = line.split('\t')
            samples = fields[9:len(fields)]
            genotype = '0'
            genotypes = []
            for s in samples:
                if s.startswith('1/1') or s.startswith('1|1') :
                    genotype = '2'
                elif s.startswith('0/1') or s.startswith('1/0') or s.startswith('0|1') or s.startswith('1|0'):
                    genotype = '1'
                elif s.startswith('0/0') or s.startswith('0|0'):
                    genotype = '0'
                else:
                    genotype = UNKNOWN

                genotypes.append(genotype)

            newline = '\t'.join(genotypes)
            fh_tmp_out.write(newline.strip()+'\n')

    fh.close()
    fh_tmp_out.close()
    print ("Tmp ped file is written to " + tmp_ped)
    print ("Transposing tmp ped file to " + ped)
    t(filein=tmp_ped, fileout=ped)
    delete (filename=tmp_ped)





###############################################################################
######################## Method to parse INFO in VCF file #####################

""" Parse key-values pairs"""
def parse_field(text, key, sep1, sep2):
    fields = text.strip().split(sep1)
    onefield=set([])
    for f in fields:
        pairs = f.split(sep2)
        if str(pairs[0]) == str(key):
        #if str(pairs[0]).find(str(key)) > -1:
            onefield.add(str(pairs[1]))

    if (len(onefield) > 0):
        return ";".join(dedup(onefield))
    else:
        return '.'
""" Parse INFO field """
def parse_info(filename, include_samples, output_file_base, extout, comment_char='##'):
    include_samples=str2bool(sys.argv[2])
    outfile = output_file_base+'.' +extout
    fh_out = open(outfile, "w")

    print "Preparing Variant file "
    print "Input file " + filename

    fh=getFh(filename)
    ok=True
    for line in fh:
        line = line.strip()
        if line.startswith(comment_char):
            #fh_out.write(line +'\n')
            print line

        elif line.startswith('CHROM') or line.startswith('#CHROM') :
            fields = line.split('\t')
            samples = []
            sample_names =[]
            if len(fields)<9:
              print ('')
              print "######## ERROR PARSING FILE ##############"
              print "## VCF file format must have at least 9 fields: "
              print "## CHROM POS     ID        REF ALT    QUAL FILTER INFO   FORMAT"
              print "## followed by sample names"
              print "## See http://www.1000genomes.org/wiki/Analysis/vcf4.0 for documentation"
              print "## You may use our parse_info_table.py script to parse INFO in tabular file format "
              print "##########################################"
              ok=False
              break

            if len(fields)>=10:
              samples = fields[9:len(fields)]
              sample_names = '\t'.join(samples).strip()


            header_head='CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER'.strip()
            header_body='\t'.join(KEYS).strip()
            header_tail='INFO\tFORMAT'.strip()

            if include_samples==True and len(fields) >= 10:
                header_tail=header_tail+ '\t' + sample_names

            header=(header_head+'\t'+header_body+'\t'+header_tail).strip()
            fh_out.write(header.strip()+'\n')

        else:
            fields = line.split('\t')
            if len(fields)>=10:
              samples = fields[9:len(fields)]
              sample_names = '\t'.join(samples).strip()

            chrom = str(fields[0])
            pos = str(fields[1])
            id = str(fields[2])
            ref = str(fields[3])
            alt  = str(fields[4])
            qual = str(fields[5])
            filter = str(fields[6])
            info = str(fields[7])
            format=str(fields[8])

            parsed_names=[]
            for k in KEYS:
                parsed_names.append(str(parse_field(info, str(k),SEP1,SEP2)))

            newline_head=chrom+'\t'+pos+'\t'+id+'\t'+ref+'\t'+alt+'\t'+qual+'\t'+filter
            newline_body='\t'.join(parsed_names).strip()
            newline_tail=(info+'\t' +format).strip()
            if include_samples == True and len(fields) >= 10:
                newline_tail=newline_tail+ '\t' + sample_names
            newline = newline_head + '\t' + newline_body + '\t' + newline_tail

            fh_out.write(newline.strip()+'\n')


    fh.close()
    fh_out.close()
    if ok:
        print ('File was saved as ' + outfile)
        print ('Deleting temporary file ' + filename)
        #delete(filename=filename)


###############################################################################
############# Start/end position for each gene     ############################

"""find the keys as a list given a value"""
def get_key(dic, value):
    return map (int, [item[0] for item in dic.items() if item[1] == value])

""" make a list of genes in the VCF file and provide index for all of them"""
def make_list_of_genes(vcffile, output_file_base, genecol, extout=GENE_EXT, geneofinterest='None', comment_char='##'):
    outfile = output_file_base+'.'+extout
    fh_out = open(outfile, "w")
    fh=getFh(filename=vcffile)
    print "Preparing Gene file "
    print "Input file " + vcffile
    gene_col_ind=-1
    eff_col_ind=-1;
    impact_col_ind=-1;
    count=1
    ugenes=[]
    genes_dict = {}
    for line in fh:
        line = line.strip()
        if line.startswith(comment_char):
            pass

        elif line.startswith('CHROM') or line.startswith('#CHROM') :
            fields = line.split('\t')
            gene_col_ind=find_first_index(fields, genecol)
            eff_col_ind = find_first_index(fields, 'SNPEFF_EFFECT')
            impact_col_ind = find_first_index(fields, 'SNPEFF_IMPACT')

            if gene_col_ind==-1:
                print genecol + ' is not found in the header of the file'
                sys.exit(2)
        else:
            fields = line.split('\t')
            gene=fields[gene_col_ind]
            eff=fields[eff_col_ind]
            imp=fields[impact_col_ind]
            #print eff
            if gene != '.' and eff != 'UPSTREAM' and eff != 'DOWNSTREAM' and imp !='MODIFIER':
                gns=gene.split(';')
                if len(gns)<2:
                    genes_dict[count]=gene
                    ugenes.append(gene)
                else:
                    genes_dict[count]=gns[0]
                    ugenes.append(gns[0])
            count=count+1

    print('')
    print ('Total genes:' + str(len(ugenes)))
    ugenes = list(set(ugenes))
    print ('Unique genes: ' + str(len(ugenes)))

    for ug in ugenes:
        keys=get_key(dic=genes_dict, value=ug)
        l= ug + '\t' + str(min(keys)) + '\t' + str(max(keys))
        fh_out.write(l.strip()+'\n')
        print(l)
    fh.close()
    fh_out.close()
    print('Result is written to ' + outfile)


###############################################################################
###############################################################################


def usage():
    print "Usage:  %s" % os.path.basename(sys.argv[0] + ' <vcf_file> <ped_file> <output_file_base>  <True\False>  <True\False>')
    print "For help:  %s" % os.path.basename(sys.argv[0] + ' -h or --help')

def main():
    if len(sys.argv)>5:
        args = sys.argv[1:]
        if "-h" in args or "--help" in args:
            usage()
            print_help()
            sys.exit(2)
        else:
            filename = sys.argv[1]
            pedfile = sys.argv[2]
            output_file_base = sys.argv[3]
            passonly = str2bool(sys.argv[4])
            changePheno = str2bool(sys.argv[5])


            for arg in sys.argv:
                print str(arg)

            if  os.path.exists(filename) and os.path.exists(pedfile) and (filename.endswith(".vcf") or filename.endswith(".vcf.gz")):
                t0 = time()

                vcf_file=filename+'.trios.tmp.vcf'
                reshufle_samples(filename=filename, pedfile=pedfile, passonly=passonly, comment_char='##')
                vcf2pedSimple(vcf_file=vcf_file, pedfile=pedfile, output_file_base=output_file_base, changePheno=changePheno, header_line_start='#CHROM', comment_chars='##', sep='\t')
                print "Applying simple conversion 1/1=2, 1/0=1, 0/0=0"
                print "All ./. are set to " + UNKNOWN


                print('')
                if os.path.exists(vcf_file) and os.path.isfile(vcf_file):
                    parse_info(filename=vcf_file, include_samples='N', output_file_base=output_file_base, extout=VAR_EXT, comment_char='##')
                    print 'File '  + vcf_file + ' does exists '
                else:
                    print 'File '  + vcf_file + ' does not exist '
                    sys.exit(1)


                ## gene file
                print('')
                vcffile=output_file_base+'.'+VAR_EXT
                make_list_of_genes(vcffile=vcffile, output_file_base=output_file_base, genecol=GENECOL)


                print('Finish')
                t1 = time()
                print 'Excecution took %f ' %(t1 - t0) + ' sec'

            else:
                print "Input files either do not exist or corrupted or do not have proper extension"
                print "Please make sure that VCF file has .vcf or .vcf.gz extension"


    elif len(sys.argv)==2:
        args = sys.argv[1:]
        if "-h" in args or "--help" in args:
            usage()
            print_help()
            sys.exit(2)
        else:
            usage()

    else:
        usage()

if __name__ == "__main__":
    sys.exit(main())
