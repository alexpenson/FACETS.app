#!/opt/common/CentOS_6-dev/python/python-2.7.10/bin/python

import argparse, os, sys, re, subprocess, itertools, errno, csv
## import cmo

SDIR = os.path.dirname(os.path.realpath(__file__))


# def make_sure_path_exists(path):
#     """
#     http://stackoverflow.com/questions/273192/in-python-check-if-a-directory-exists-and-create-it-if-necessary
#     """
#     try:
#         os.makedirs(path)
#     except OSError as exception:
#         if exception.errno != errno.EEXIST:
#             raise

def make_sure_path_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory) 

        
def slugify(value):
    """
    Normalizes string, removes non-alpha characters,
    and converts spaces to hyphens.
    http://stackoverflow.com/questions/295135/turn-a-string-into-a-valid-filename-in-python
    """
    ##    import unicodedata
    ##    value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore')
    value = re.sub('[^\w\s-]', '', value).strip()
    re.sub('[-\s]+', '-', value)
    return(value)


def fromcounts(args):
    print "run " + SDIR + "/bin/doFacets.R"

def run(args):
    """read pairs file and for each line do:
    1) counts SNP position in tumor and normal bam files
    2) merge two counts files with threshold on normal depth
    3) run facets on merged counts"""
    
    ### check for Matched_Norm_Sample_Barcode, if column exists then use it...
        
    cmd_list = list()
    ### LOOP OVER ENTRIES IN PAIRS FILE
    pairs_dict = csv.DictReader(args.pairs_file, delimiter='\t')
    for facets_run in pairs_dict:
        Tumor_Sample_Barcode = slugify(facets_run['Tumor_Sample_Barcode'])
        t_bamfile = facets_run['t_bamfile']
        n_bamfile = facets_run['n_bamfile']
        # print Tumor_Sample_Barcode
        # print t_bamfile
        # print n_bamfile
        # print SDIR

        ### CHECK FOR EXISTENCE OF INPUT FILES

        ### GENERATE NAMES FOR COUNTS FILES
        n_countsfile = "counts/n_counts____" + Tumor_Sample_Barcode + ".dat.gz"
        t_countsfile = "counts/t_counts____" + Tumor_Sample_Barcode + ".dat.gz"
        countsMerged_file = "counts/countsMerged____" + Tumor_Sample_Barcode + ".dat.gz"


        ### CHECK FOR EXISTENCE OF COUNTS FILES
        # with gzip.open('file.txt.gz', 'rb') as f:
        #     file_content = f.read()


        ### BUILD LIST OF COMMANDS TO BE EXECUTED (WHILE WAITING FOR ERRORS)
                
        ### GET BASE COUNTS
        
        counts_cmd = ('bsub -We 59 -o Log/ -e Err/ -J %s_count_%s -R "rusage[mem=40]" ' 
                      '%s/getBaseCountsZZAutoWithName.sh %s %s')
        
        n_counts_cmd = counts_cmd % ("n", Tumor_Sample_Barcode, SDIR, n_countsfile, n_bamfile)
        cmd_list.append(n_counts_cmd)
        t_counts_cmd = counts_cmd % ("t", Tumor_Sample_Barcode, SDIR, t_countsfile, t_bamfile)        
        cmd_list.append(t_counts_cmd)
        
        
        ### MERGE COUNTS
        wait_string = ''
        ##        wait_string = '-w "post_done(*_count_%s)"' % (Tumor_Sample_Barcode)
        
        merge_cmd = ('bsub -We 59 -o Log/ -e Err/ -n 2 -R "rusage[mem=60]" -J merge_%s %s '
                     '"%s/mergeTN.R %s %s | gzip -9 -c > %s"')
        merge_cmd = merge_cmd % (Tumor_Sample_Barcode, wait_string, SDIR, t_countsfile, n_countsfile, countsMerged_file)
        cmd_list.append(merge_cmd)

        
        ### EXECUTE COMMANDS

        ### make output directories for the project
        make_sure_path_exists("counts")
        make_sure_path_exists("cval_50")

        for facets_run in pairs_dict:
            ### MAKE OUTPUT FOLDER FOR EACH facets_run
            make_sure_path_exists("cval_50/" + Tumor_Sample_Barcode)
            
        for cmd in cmd_list:
            # ###            subprocess.call(cmd) ### why doesn't this work?
            subprocess.call(["echo", cmd])
            #os.system(cmd)


def call(args):
    print "call genes"

def merge(args):
    print "run " + SDIR + "/bin/postFacets.sh"
    
def maf(args):
    print "run JoinFACETS2maf"


        

        

if __name__ =='__main__':

    ### ARGUMENTS

    parser = argparse.ArgumentParser(description="run FACETS analysis")
    subparsers = parser.add_subparsers(help='sub-command help')

    ### facets.py run
    parser_run = subparsers.add_parser('run', help='run FACETS from bam files')
    parser_run.add_argument('pairs_file', type=argparse.FileType('r'), 
                            help=('Tumor/Normal pairs file: must contain columns '  
                                  'Tumor_Sample_Barcode, t_bamfile & n_bamfile, tab-delimited'))
    parser_run.set_defaults(func=run)
    ### optional aruments: cval ndepth etc.    

    ### facets.py fromcounts
    parser_fromcounts = subparsers.add_parser('fromcounts', help='extract SNP fromcounts from bam files')
    parser_fromcounts.add_argument('counts_files', nargs = '*', type=argparse.FileType('r'), 
                               help=('Tumor/Normal pairs file: must contain columns '  
                                     'Tumor_Sample_Barcode, t_bamfile & n_bamfile, tab-delimited'))
    parser_fromcounts.set_defaults(func=fromcounts)
    ### optional aruments: cval ndepth etc.    

    args = parser.parse_args()
    args.func(args)
    

