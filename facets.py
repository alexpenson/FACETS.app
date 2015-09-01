#!/opt/common/CentOS_6-dev/python/python-2.7.10/bin/python

import argparse, os, sys, re, subprocess, itertools, errno, csv, gzip
## import cmo

SDIR = os.path.dirname(os.path.realpath(__file__))


def make_sure_path_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory) 

def gzip_file_with_size(file_path):
    if not os.path.exists(file_path): return(False)
    with gzip.open(file_path, 'rb') as f:
        for i, l in enumerate(f):
            i = i+1
            if i > 10: return(True)
        return(False)

        
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


def runlsf(args):
    """read pairs file and for each line create LSF commands to:
    1) counts SNP position in tumor and normal bam files
    2) merge two counts files with threshold on normal depth
    3) run facets on merged counts"""

    print args
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
        n_countsfile_exists = gzip_file_with_size(n_countsfile)
        t_countsfile_exists = gzip_file_with_size(t_countsfile)
        countsMerged_file_exists = gzip_file_with_size(countsMerged_file)


        ### BUILD LIST OF COMMANDS TO BE EXECUTED (WHILE WAITING FOR ERRORS)
                
        ### GET BASE COUNTS

        wait_string = '' ### string to add to later bsub commands

        counts_cmd = ('bsub -We 59 -o LSF/ -e Err/ -J %s_count_%s -R "rusage[mem=40]" ' 
                      '%s/bin/getBaseCountsZZAutoWithName.sh %s %s')
        
        if not n_countsfile_exists and not countsMerged_file_exists:
            n_counts_cmd = counts_cmd % ("n", Tumor_Sample_Barcode, SDIR, n_countsfile, n_bamfile)
            cmd_list.append(n_counts_cmd)
            wait_string = '-w "post_done(n_count_%s) && post_done(t_count_%s)"' % (Tumor_Sample_Barcode, Tumor_Sample_Barcode)
            
        if not t_countsfile_exists and not countsMerged_file_exists:
            t_counts_cmd = counts_cmd % ("t", Tumor_Sample_Barcode, SDIR, t_countsfile, t_bamfile)        
            cmd_list.append(t_counts_cmd)
            wait_string = '-w "post_done(n_count_%s) && post_done(t_count_%s)"' % (Tumor_Sample_Barcode, Tumor_Sample_Barcode)
                    
        
        ### MERGE COUNTS
        if not countsMerged_file_exists:
            merge_cmd = ('bsub -We 59 -o LSF/ -e Err/ -n 2 -R "rusage[mem=60]" -J merge_%s %s '
                         '"%s/bin/mergeTN.R %s %s | gzip -9 -c > %s"')
            merge_cmd = merge_cmd % (Tumor_Sample_Barcode, wait_string, SDIR, t_countsfile, n_countsfile, countsMerged_file)
            cmd_list.append(merge_cmd)
            wait_string = '-w "post_done(merge_%s)"' % (Tumor_Sample_Barcode)
        
        # ### NORM COUNTS
        # norm_cmd = ('bsub -We 59 -o LSF/ -e Err/ -n 2 -R "rusage[mem=60]" -J norm_%s %s '
        #              '"%s/bin/norm_normal_depth.R %s | gzip -9 -c > %s"')
        # norm_cmd = norm_cmd % (Tumor_Sample_Barcode, wait_string, SDIR, countsMerged_file, countsMerged_file_norm)
        # cmd_list.append(norm_cmd)


        ### DO FACETS
        doFacets_cmd = ('bsub -We 59 -o LSF/ -e Err/ -J facets_%s %s '
                        '%s/bin/doFacets.R %s %s %s %s')
        doFacets_cmd = doFacets_cmd % (Tumor_Sample_Barcode,
                                       wait_string,
                                       SDIR,
                                       " ".join(args.facets_args),
                                       countsMerged_file,
                                       Tumor_Sample_Barcode,
                                       args.outputdir)
        cmd_list.append(doFacets_cmd)
        
    ### EXECUTE COMMANDS

    ### make output directories for the project
    make_sure_path_exists("counts")
    make_sure_path_exists(args.outputdir)

    for facets_run in pairs_dict:
        ### MAKE OUTPUT FOLDER FOR EACH facets_run
        make_sure_path_exists("cval_50/" + Tumor_Sample_Barcode)
        
    for cmd in cmd_list:
        subprocess.call(["echo", cmd])
        subprocess.call(cmd, shell=True)
        #os.system(cmd)

def fromcounts(args):
    ### BUILD COMMAND LIST
    for countsMerged_file in args.counts_files:
        doFacets_cmd = ('bsub -We 59 -o LSF/ -e Err/ -J facets_%s %s '
                        '%s/bin/doFacets.R %s %s')
        doFacets_cmd = doFacets_cmd % (Tumor_Sample_Barcode, wait_string, SDIR, " ".join(args.facets_args), countsMerged_file)
        cmd_list.append(doFacets_cmd)

    ### EXECUTE COMMAND LIST
    for cmd in cmd_list:
        subprocess.call(["echo", cmd])
        subprocess.call(cmd, shell=True)

# def norm(args):
#     print "run " + SDIR + "/bin/norm_normal_depth.R"
            
def check(args):
    print "check output files"

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

    ### ./facets.py run
    parser_runlsf = subparsers.add_parser('runlsf', help='create LSF commands to run FACETS from bam files')
    parser_runlsf.add_argument('-o', '--outputdir',action='store', default='output',
                               help='directory for output files')
    parser_runlsf.add_argument('pairs_file', type=argparse.FileType('r'), 
                               help=('Tumor/Normal pairs file: must contain columns '  
                                     'Tumor_Sample_Barcode, t_bamfile & n_bamfile, tab-delimited'))
    parser_runlsf.add_argument('facets_args', nargs=argparse.REMAINDER,
                               help='remaining arguments are sent to doFacets.R') 
    parser_runlsf.set_defaults(func=runlsf)

    
    ### ./facets.py fromcounts
    parser_fromcounts = subparsers.add_parser('fromcounts', help='run FACETS from merged counts files')
    parser_fromcounts.add_argument('counts_files', nargs = '*', type=argparse.FileType('r'), 
                               help='merge counts files counts/countsMerged__*.dat.gz')
    parser_fromcounts.add_argument('facets_args', nargs=argparse.REMAINDER,
                                   help='remaining arguments are sent to doFacets.R') 
    parser_fromcounts.set_defaults(func=fromcounts)

    # ### ./facets.py norm
    # parser_norm = subparsers.add_parser('norm', help='run FACETS from merged counts files')
    # parser_norm.add_argument('counts_files', nargs = '*', type=argparse.FileType('r'), 
    #                            help='merge counts files counts/countsMerged__*.dat.gz')
    # ### remaining arguments are sent to doFacets.R
    # parser_norm.set_defaults(func=norm)


    args = parser.parse_args()
##    print args
    args.func(args)
    

