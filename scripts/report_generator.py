import getopt, sys, os
from optparse import OptionParser
import os
import subprocess
import ntpath
import zipfile
from bdbag import bdbag_api
from shutil import copy

def main():
    parser = OptionParser()
    parser.add_option("-b", "--bdbag", dest="bdbag",help="Path to bdbag archive", metavar="BDBAG")
    parser.add_option("-o", "--outdir", dest="outdir",help="Path to Output directory", metavar="PATH")
    parser.add_option("-n", "--name", dest="pname",help="Name of Project must match that when generating bag", metavar="NAME")
    parser.add_option("-1", "--fastqc", dest="fqc",help="Make fastQC Report", metavar="FQC", action='store_true')
    parser.add_option("-2", "--align", dest="aln",help="Make Alignment Report", metavar="ALN", action='store_true')
    parser.add_option("-3", "--ge", dest="ge",help="Make GE Report", metavar="GE", action='store_true')
    parser.add_option("-4", "--de", dest="de",help="Make DE Report", metavar="DE", action='store_true')
    parser.add_option("-a", "--all", dest="all",help="Make All Reports", metavar="ALL", action='store_true')
    parser.add_option("-p", "--prok", dest="prok",help="Make Prok Reports", metavar="PROK", action='store_true')
    parser.add_option("-u", "--update", dest="update",help="Update BdBag", metavar="Update", action='store_true')
    parser.add_option("-i", "--info", dest="info",help="Path to Info file", metavar="PATH")
    parser.add_option("-m", "--mapping", dest="mapping",help="Path to mapping file", metavar="PATH")

    (options, args) = parser.parse_args()
    #wrap_FQC="/usr/local/packages/report_generation/wrapper_FastQC.R"
    
    wrap_FQC="/home/apaala.chatterjee/RNA_Report/wrapper_FastQC.R"
    wrap_ALN_prok="/home/apaala.chatterjee/RNA_Report/wrapper_Alignment_prok.R"
    wrap_ALN= "/home/apaala.chatterjee/RNA_Report/wrapper_Alignment.R"
    wrap_GE="/home/apaala.chatterjee/RNA_Report/wrapper_GE.R"
    wrap_DE="/home/apaala.chatterjee/RNA_Report/wrapper_DE.R"
    wrap_DE_mapping="/home/apaala.chatterjee/RNA_Report/wrapper_DE_mapping.R"
    wrap_GE_mapping="/home/apaala.chatterjee/RNA_Report/wrapper_GE_mapping.R"
    counts_script="/usr/local/packages/report_generation/Generate_all_counts.R"
    rpath="/usr/local/packages/r-3.4.0/bin/Rscript"

    if options.bdbag and options.outdir and options.pname:
        out_base = os.path.normpath(options.outdir)
        bdbag_api.extract_bag(options.bdbag, output_path = out_base)
        extracted_path = os.path.normpath(out_base + "/" + options.pname + "/data/")
        copy(options.info, extracted_path)

        if(options.fqc or options.all):
            generate_fastqc_report(extracted_path, wrap_FQC, rpath, options.pname, options.info)

        if((options.aln or options.all) and not options.prok):
            generate_alignment_report(extracted_path, wrap_ALN, rpath, options.pname, options.info)
    
        if((options.aln or options.all) and options.prok):
            generate_alignment_report_prok(extracted_path, wrap_ALN_prok, rpath, options.pname, options.info)
    
        if((options.ge or options.all) and not options.mapping):
            generate_ge_report(extracted_path, wrap_GE, rpath, options.pname, options.info)

        elif((options.ge or options.all) and options.mapping):
            print("In mapping file GE")
            map_to_GE(extracted_path, wrap_GE_mapping, rpath, options.pname, options.info, options.mapping)
        
        if((options.de or options.all) and not options.mapping):
            generate_de_report(extracted_path, wrap_DE, rpath, options.pname, options.info)

        elif((options.de or options.all) and options.mapping):
            print("In mapping file DE")
            map_to_DE(extracted_path, wrap_DE_mapping, rpath, options.pname, options.info, options.mapping)
        
        if(options.update):
            bag_path = os.path.normpath(out_base + "/" + options.pname)
            bdbag_api.make_bag(bag_path, update = True)
            bdbag_api.archive_bag(bag_path, "zip")
    
def generate_all_counts(path_to_counts):
    """
    Input: path to directory with counts files 
    
    Output: Single merged counts dataframe
    """

    counts_list = [f for f in os.listdir(path_to_counts) if f.endswith('.counts') or f.endswith('.counts.txt')]
    appended_counts_list = prepend(counts_list, path_to_counts)
    all_counts_merge = pandas.read_csv(appended_counts_list[0], sep="\t", header = None)
    count_name_0 = ntpath.basename(appended_counts_list[0]).split('.',1)[0]
    all_counts_merge.rename(columns = {1: count_name_0}, inplace = True)
    print(all_counts_merge.head())
    for i in range(1,len(appended_counts_list)):
        hold_table = pandas.read_csv(appended_counts_list[i], sep="\t", header = None)
        hold_name = ntpath.basename(appended_counts_list[i]).split('.',1)[0]
        hold_table.rename(columns = {1: hold_name}, inplace = True)
        all_counts_merge = pandas.merge(all_counts_merge, hold_table, on = 0)
    all_counts_merge.rename(columns = {0: "ID"}, inplace = True)
    print(all_counts_merge.head())
    return(all_counts_merge)

def generate_alignment_report(outdir, wrapper, rpath, project_name, info_file):
    """
    Input: path to output directory with extracted files, path to wrapper R script, path to R, name of project and path to info file. 
    
    Output: Alignment report generated and output plots saved in AlignmentFiles
    """

    output_files = os.path.normpath(outdir + "/AlignmentFiles/")
    if not os.path.exists(output_files):
        os.makedirs(output_files)
    try:
        subprocess.check_call([rpath, wrapper, project_name, outdir ,info_file, outdir], shell = False)
    except subprocess.CalledProcessError as e:
        print(e)

def generate_alignment_report_prok(outdir, wrapper, rpath, project_name, info_file):
    """
    Input: path to output directory with extracted files, path to wrapper R script, path to R, name of project and path to info file. 
    
    Output: Prok alignment report generated and output plots saved in AlignmentFiles
    """

    output_files = os.path.normpath(outdir + "/AlignmentFiles/")
    if not os.path.exists(output_files):
        os.makedirs(output_files)
    try:
        subprocess.check_call([rpath, wrapper, project_name, outdir ,info_file, outdir], shell = False)
    except subprocess.CalledProcessError as e:
        print(e)

def generate_de_report(outdir, wrapper, rpath, project_name, info_file):
    """
    Input: path to output directory with extracted files, path to wrapper R script, path to R, name of project and path to info file.

    Output: DE report generated and outputs saved in DE_Outputs
    """

    output_files = os.path.normpath(outdir + "/DE_Outputs/")
    if not os.path.exists(output_files):
        os.makedirs(output_files)
    all_count_path = os.path.normpath(outdir + "/all_counts.txt")
    try:
        subprocess.check_call([rpath, wrapper, project_name, outdir ,info_file, outdir], shell = False)
    except subprocess.CalledProcessError as e:
        print(e)
    
def generate_fastqc_report(outdir, wrapper, rpath, project_name, info_file):
    """
    Input: path to output directory with extracted files, path to wrapper R script, path to R, name of project and path to info file. 

    Output: FastQC report generate and output files FastQC_Outputs
    """

    output_files = os.path.normpath(outdir + "/FastQC_Outputs/")
    if not os.path.exists(output_files):
        os.makedirs(output_files)
    try:
        subprocess.check_call([rpath, wrapper, project_name, outdir ,info_file, outdir], shell = False)
    except subprocess.CalledProcessError as e:
        print(e)
    #syscmd=rpath +" "+ wrap_FQC +" "+ options.pname +" "+ pdir +" " + options.info + " "+pdir

def generate_ge_report(outdir, wrapper, rpath, project_name, info_file):
    """
    Input: path to output directory with extracted files, path to wrapper R script, path to R, name of project and path to info file.

    Output: GE report generate and output files in GE_Outputs 
    """

    output_files = os.path.normpath(outdir + "/GE_Outputs/")
    all_count_path = os.path.normpath(outdir + "/all_counts.txt")
    count_file_expected_path = os.path.normpath(outdir + "/Counts/")
    if not os.path.exists(output_files):
        os.makedirs(output_files)
    if not os.path.exists(all_count_path):
        merged_counts = generate_all_counts(count_file_expected_path)
        merged_counts.to_csv(path_or_buf = all_count_path, header = True, sep = "\t", index = False)
    try:
        subprocess.check_call([rpath, wrapper, project_name, outdir ,info_file, outdir], shell = False)
    except subprocess.CalledProcessError as e:
        print(e)

def map_to_DE(outdir, wrapper, rpath, project_name, info_file, mappingf):
    """
    Input: path to output directory with extracted files, path to wrapper R script, path to R, name of project and path to info file.

    Output: DE report generate and output files in DE_Outputs 
    """

    output_files = os.path.normpath(outdir + "/DE_Outputs/")
    if not os.path.exists(output_files):
        os.makedirs(output_files)
    all_count_path = os.path.normpath(outdir + "/all_counts.txt")
    try:
        subprocess.check_call([rpath, wrapper, project_name, outdir ,info_file, outdir, mappingf], shell = False)
    except subprocess.CalledProcessError as e:
        print(e)

def map_to_GE(outdir, wrapper, rpath, project_name, info_file, mappingf):
    """
    Input: path to output directory with extracted files, path to wrapper R script, path to R, name of project, path to info file and path to mapping
           file

    Output: GE report generate and output files in GE_Outputs 
    """

    output_files = os.path.normpath(outdir + "/GE_Outputs/")
    if not os.path.exists(output_files):
        os.makedirs(output_files)
    all_count_path = os.path.normpath(outdir + "/all_counts.txt")
    try:
        subprocess.check_call([rpath, wrapper, project_name, outdir ,info_file, outdir, mappingf], shell = False)
    except subprocess.CalledProcessError as e:
        print(e)

if __name__ == '__main__':
    main()
