#!/opt/common/CentOS_6-dev/python/python-2.7.10/bin/python

import argparse, os, sys, re, subprocess, itertools, errno, csv
## import cmo

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


if __name__ =='__main__':

    ### ARGUMENTS

    parser = argparse.ArgumentParser(description="run FACETS")
    parser.add_argument('pairs_file', type=argparse.FileType('r'), 
                        help='''Tumor/Normal pairs file: must contain columns Tumor_Sample_Barcode, t_bamfile & n_bamfile, tab-delimited''')
    ### optional aruments: cval ndepth etc.    

    args = parser.parse_args()

    ### check for Matched_Norm_Sample_Barcode, if column exists then use it...


    
    ### OTHER VARIABLES

    SDIR = os.path.dirname(os.path.realpath(__file__))

    

    ### DO STUFF
    
    ### make output directories for the project
    make_sure_path_exists("counts")
    make_sure_path_exists("cval_50")


    cmd_list = list()
    ### LOOP OVER ENTRIES IN PAIRS FILE
    pairs_dict = csv.DictReader(args.pairs_file, delimiter='\t')
    for facets_run in pairs_dict:
        Tumor_Sample_Barcode = slugify(facets_run['Tumor_Sample_Barcode'])
        t_bamfile = facets_run['t_bamfile']
        n_bamfile = facets_run['n_bamfile']
        print Tumor_Sample_Barcode
        print t_bamfile
        print n_bamfile
        print SDIR

        ### MAKE OUTPUT FOLDER FOR EACH facets_run

        make_sure_path_exists("cval_50/" + Tumor_Sample_Barcode)
        
        ### GENERATE NAMES FOR OUTPUT FILES
        n_countsfile = "counts/n_counts____" + Tumor_Sample_Barcode + ".dat.gz"
        t_countsfile = "counts/t_counts____" + Tumor_Sample_Barcode + ".dat.gz"
        countsMerged_file = "counts/countsMerged____" + Tumor_Sample_Barcode + ".dat.gz"


        ### CHECK FOR EXISTENCE OF OUTPUT FILES
        # with gzip.open('file.txt.gz', 'rb') as f:
        #     file_content = f.read()

        
        ### GET BASE COUNTS
        
        counts_cmd = 'bsub -We 59 -o Log/ -e Err/ -J %s_count_%s -R "rusage[mem=40]" %s/getBaseCountsZZAutoWithName.sh %s %s'

        n_counts_cmd = counts_cmd % ("n", Tumor_Sample_Barcode, SDIR, n_countsfile, n_bamfile)
        t_counts_cmd = counts_cmd % ("t", Tumor_Sample_Barcode, SDIR, t_countsfile, t_bamfile)

        print n_counts_cmd
        cmd_list.append(n_counts_cmd)
        print t_counts_cmd
        cmd_list.append(t_counts_cmd)


        ### MERGE COUNTS
        wait_string = ''
##        wait_string = '-w "post_done(*_count_%s)"' % (Tumor_Sample_Barcode)

        merge_cmd = '''
        bsub -We 59 -o Log/ -e Err/ -n 2 -R "rusage[mem=60]" -J merge_%s %s \
        "%s/mergeTN.R %s %s | gzip -9 -c > %s"
        ''' % (Tumor_Sample_Barcode, wait_string, SDIR, t_countsfile, n_countsfile, countsMerged_file)

        print merge_cmd
        cmd_list.append(merge_cmd)

        ### build list of commands to be executed (while waiting for errors)

        
        ### execute commands
        for cmd in cmd_list:
###            subprocess.call(cmd) ### why doesn't this work?
            os.system(cmd)
