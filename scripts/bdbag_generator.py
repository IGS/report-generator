import getopt, sys, os
from optparse import OptionParser
import os
import subprocess
import ntpath
import zipfile
import argparse
import shutil
import pathlib
import pandas
from bdbag import bdbag_api
import bagit

def main():
    parser = argparse.ArgumentParser( description='BDBag Generator for Report Generation Tool')
    parser.add_argument("-d", "--dir", required=True, dest="pdir",help="Path to Ergatis Structure directory", metavar="PATH")
    parser.add_argument("-o", "--outdir", required=True, dest="outdir",help="Path to Output directory", metavar="PATH")
    parser.add_argument("-n", "--name", required=True, dest="pname",help="Name of Project must match that when generating bag", metavar="NAME")
    parser.add_argument("-1", "--fastqc", dest="fqc",help="Make fastQC Report", action='store_true')
    parser.add_argument("-2", "--align", dest="aln",help="Make Alignment Report", action='store_true')
    parser.add_argument("-3", "--ge", dest="ge",help="Make GE Report", action='store_true')
    parser.add_argument("-4", "--de", dest="de",help="Make DE Report", action='store_true')
    parser.add_argument("-a", "--all", dest="allreports",help="Make All Reports", action='store_true')
    parser.add_argument("-u", "--update", dest="update",help="Update BdBag", action='store_true')
    parser.add_argument("-p", "--pid", required=True, dest="pid",help="Pipeline ID")
    parser.add_argument("-c", "--cloud", dest="cloud",help="Pull files for Cloud")
    parser.add_argument("-b", "--bam", dest="bam",help="Pull BAMs into BAG", action='store_true')
    parser.add_argument("-s", "--sam", dest="sam",help="Pull SAMs into BAG", action='store_true')
    parser.add_argument("-m", "--make", dest="make",help="Make BDBAG", action='store_true')
    args = parser.parse_args()
    
    ergatis_repository = os.path.normpath(args.pdir + "/output_repository/")

    #Make output directory with name
    output_dir = os.path.normpath(args.outdir + "/" + args.pname)

    #Check What reports are requested and pulling relevant files
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if args.fqc and not args.allreports:
        print("Starting FastQC file gathering")
        only_fqc = generate_fastqc_report(output_dir, ergatis_repository, args.pid)
    if args.aln and not args.allreports:
        print("Starting Alignment file gathering")
        only_aln = generate_alignment_report(output_dir, ergatis_repository, args.pid)
    if args.ge and not args.allreports:
        print("Starting GE file gathering")
        only_ge = generate_ge_report(output_dir, ergatis_repository, args.pid)
    if args.de and not args.allreports:
        only_de = generate_de_report(output_dir, ergatis_repository, args.pid)
    if args.allreports:
        generate_all_reports(output_dir, ergatis_repository, args.pid)
    if args.make and not args.update:
        create_bag(output_dir, False)
    if args.update:
        create_bag(output_dir, True)

def copy_bam_files():
    pass

def copy_cloud_files(outdir, ergatis_base, ergatis_pid):
   pass 

def copy_files_to_dir(these_files, to_here):
    for i in range(len(these_files)):
        syscmd = "cp "+ these_files[i] + " " + to_here[i]
        os.system(syscmd) 

def create_bag(output_dir, update):
    """Create/Update and archive a BDBag from the contents of a passed-in directory."""
    bdbag_api.make_bag(output_dir, update=update)
    return bdbag_api.archive_bag(output_dir, "zip")

def generate_all_counts(path_to_counts):
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

def generate_all_reports(outdir, ergatis_repository, ergatis_pid):
    """Generate all possible reports."""
    only_fqc = generate_fastqc_report(output_dir, ergatis_repository, ergatis_pid)
    only_aln = generate_alignment_report(output_dir, ergatis_repository, ergatis_pid)
    only_ge = generate_ge_report(output_dir, ergatis_repository, ergatis_pid)
    only_de = generate_de_report(output_dir, ergatis_repository, ergatis_pid)

def generate_alignment_report(outdir, ergatis_repository, ergatis_pid):
    """
    Input: path to input and output directory 
    
    Output: Alignment files copied
    """
    print("Copying Alignment Files")
    aln_path = os.path.normpath(ergatis_repository + "/wrapper_align/" + ergatis_pid +"_wrap/Summary.txt")
    syscmd = "cp " + aln_path + " " + outdir
    print(syscmd)
    os.system(syscmd)

def generate_de_report(outdir, ergatis_repository, ergatis_pid):
    de_paths = [os.path.normpath(ergatis_repository + "/deseq/" + ergatis_pid + "_differential_expression/i1/g*/*.counts.txt"), os.path.normpath(ergatis_repository + "/deseq/"+ ergatis_pid + "_differential_expression/i1/g*/*de_genes.txt")]
    de_dir = os.path.normpath(outdir + "/DE/")
    if not os.path.exists(de_dir):
        os.makedirs(de_dir)
    copy_files_to_dir(de_paths, [de_dir,de_dir])
    counts_default = os.path.normpath(outdir + "/" +"all_counts.txt")
    counts_de = [f for f in os.listdir(de_dir) if f.endswith('.counts.txt')]
    count_file_de = os.path.normpath(de_dir + "/" + counts_de[0])
    print(count_file_de)
    if not os.path.exists(counts_default):
        print("Copying all counts")
        shutil.copy(count_file_de , counts_default)  


def generate_ge_report(outdir, ergatis_repository, ergatis_pid):
    """
    Input: path to input and output directory 
    
    Output: GE files copied
    """
    print("In GE")
    ge_paths = [os.path.normpath(ergatis_repository + "/rpkm_coverage_stats/" + ergatis_pid + "_rpkm_cvg/i1/g*/genic_coverage/*.txt") , os.path.normpath(ergatis_repository + "/htseq/" + ergatis_pid + "*_counts/i1/g*/*.counts")]    
    print(ge_paths)
    ge_dir = ["/RPKM/","/Counts/"]
    ge_full_paths = prepend(ge_dir, outdir)
    print(ge_full_paths)
    makedir(ge_full_paths)
    copy_files_to_dir(ge_paths, ge_full_paths)
    all_counts = generate_all_counts( ge_full_paths[1])    
    counts_out = os.path.normpath(outdir + "/" +"all_counts.txt")
    all_counts.to_csv(path_or_buf = counts_out, header = True, sep = "\t", index = False)

def generate_fastqc_report(outdir, ergatis_repository, ergatis_pid):
    """
    Input: path to input and output directory 
    
    Output: Path and Name to FastQC report
    """
    fqc_path = os.path.normpath(ergatis_repository + "/fastqc_stats/" + ergatis_pid + "_fastqc/i1/g*/*/Images/")
    print("Staring Directory Setup")

    dir_list = ["/FastQC_Files", "/FastQC_Files/KmerProfiles", "/FastQC_Files/BaseQuality", "/FastQC_Files/adapter_content", "/FastQC_Files/duplication_levels", "/FastQC_Files/per_base_n_content", "/FastQC_Files/per_base_sequence_content", "/FastQC_Files/per_sequence_gc_content", "/FastQC_Files/per_sequence_quality", "/FastQC_Files/per_tile_quality", "/FastQC_Files/sequence_length_distribution"] 

    extension_list = ["/*_sequence.kmer_profiles.png", "/*_base_quality.png", "/*_sequence.adapter_content.png", "/*_sequence.duplication_levels.png", "/*_sequence.per_base_n_content.png", "/*_sequence.per_base_sequence_content.png", "/*_sequence.per_sequence_gc_content.png", "/*_sequence.per_sequence_quality.png", "/*_sequence.per_tile_quality.png", "/*_sequence.sequence_length_distribution.png"]

    copy_to = ["/FastQC_Files/KmerProfiles", "/FastQC_Files/BaseQuality", "/FastQC_Files/adapter_content", "/FastQC_Files/duplication_levels", "/FastQC_Files/per_base_n_content", "/FastQC_Files/per_base_sequence_content", "/FastQC_Files/per_sequence_gc_content", "/FastQC_Files/per_sequence_quality", "/FastQC_Files/per_tile_quality", "/FastQC_Files/sequence_length_distribution"]

    appended_dir_list = prepend(dir_list, outdir)
    appended_copy_to = prepend(copy_to, outdir)
    appended_to_be_copied = prepend(extension_list, fqc_path)
    makedir(appended_dir_list) 
    copy_files_to_dir(appended_to_be_copied, appended_copy_to)
        

def makedir(pathlist):
    """
    Input: List of FULL Paths to be made if they dont exist 
    
    Output: Confirmation if the paths were made
    """
    for path in pathlist:
        if not os.path.exists(path):
            os.makedirs(path)
    return ""
    
def prepend(list, str):
    """
    Input: List of files in the input directory and path to input directory

    Output: List of full paths to log files in input directory
    """
    str += '{0}'
    list = [str.format(i) for i in list]
    return(list)

if __name__ == '__main__':
    main()
