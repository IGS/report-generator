import sys, os
import argparse
import shutil
import pandas
from bdbag import bdbag_api
import bagit
import glob

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
    parser.add_argument("-r", "--decounts", dest="rename", help = "Path to R script to copy and rename counts files from deseq", metavar="Path")
    args = parser.parse_args()

    ergatis_repository = os.path.normpath(args.pdir + "/output_repository/")

    #Make output directory with name
    output_dir = os.path.normpath(args.outdir + "/" + args.pname)

    #Rename Counts 
    if args.rename:
        renameDE = args.rename
    else:
        renameDE= "/home/apaala.chatterjee/RNA_Report/rename_de.R"
    #print(renameDE)
    #Check What reports are requested and pulling relevant files
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if args.allreports:
        generate_all_reports(output_dir, ergatis_repository, args.pid, renameDE)
    else:
        if args.fqc:
            print("Starting FastQC file gathering")
            only_fqc = generate_fastqc_report(output_dir, ergatis_repository, args.pid)
        if args.aln:
            print("Starting Alignment file gathering")
            only_aln = generate_alignment_report(output_dir, ergatis_repository, args.pid)
        if args.ge:
            print("Starting GE file gathering")
            only_ge = generate_ge_report(output_dir, ergatis_repository, args.pid)
        if args.de:
            only_de = generate_de_report(output_dir, ergatis_repository, args.pid, renameDE)
    if args.make and not args.update:
        create_bag(output_dir, False)
    if args.update:
        create_bag(output_dir, True)

def copy_bam_files():
    pass

def copy_cloud_files(outdir, ergatis_base, ergatis_pid):
   pass

def copy_de_files(path_to_counts_basedir, DEPath, renameDE):
    ##Change to not use R script later and make it a independent function in python
    count_cmd = "Rscript "+ renameDE + " " + path_to_counts_basedir + " " + DEPath
    print(count_cmd)
    os.system(count_cmd)

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
    count_name_0 = os.path.basename(appended_counts_list[0]).split('.',1)[0]
    all_counts_merge.rename(columns = {1: count_name_0}, inplace = True)
    print(all_counts_merge.head())
    for i in range(1,len(appended_counts_list)):
        hold_table = pandas.read_csv(appended_counts_list[i], sep="\t", header = None)
        hold_name = os.path.basename(appended_counts_list[i]).split('.',1)[0]
        hold_table.rename(columns = {1: hold_name}, inplace = True)
        all_counts_merge = pandas.merge(all_counts_merge, hold_table, on = 0)
    all_counts_merge.rename(columns = {0: "ID"}, inplace = True)
    print(all_counts_merge.head())
    return(all_counts_merge)

def generate_all_reports(output_dir, ergatis_repository, ergatis_pid, renameDE):
    """Generate all possible reports."""
    only_fqc = generate_fastqc_report(output_dir, ergatis_repository, ergatis_pid)
    only_aln = generate_alignment_report(output_dir, ergatis_repository, ergatis_pid)
    only_ge = generate_ge_report(output_dir, ergatis_repository, ergatis_pid)
    only_de = generate_de_report(output_dir, ergatis_repository, ergatis_pid, renameDE)

def generate_alignment_report(outdir, ergatis_repository, ergatis_pid):
    """
    Input: path to input and output directory

    Output: Alignment files copied
    """
    # Remove redundant references and separators in path, if applicable
    outdir = os.path.normpath(outdir)
    ergatis_repository = os.path.normpath(ergatis_repository)

    print("Copying Alignment Files")
    aln_path = os.path.join(ergatis_repository, "wrapper_align", ergatis_pid +"_wrap/Summary.txt")
    syscmd = "cp " + aln_path + " " + outdir
    print(syscmd)
    os.system(syscmd)

def generate_de_report(outdir, ergatis_repository, ergatis_pid, renameDE):
    # Remove redundant references and separators in path, if applicable
    outdir = os.path.normpath(outdir)
    ergatis_repository = os.path.normpath(ergatis_repository)
    base_counts = os.path.join(ergatis_repository, "deseq", ergatis_pid + "_differential_expression/i1/g*/")
    print("In DE")
    #print(renameDE)
    de_paths = [os.path.join(ergatis_repository, "deseq", ergatis_pid + "_differential_expression/i1/g*/*.counts.txt"), os.path.join(ergatis_repository, "deseq", ergatis_pid + "_differential_expression/i1/g*/*de_genes.txt")]
    #print(de_paths)
    de_dir = os.path.join(outdir, "DE")
    if not os.path.isdir(de_dir):
        os.makedirs(de_dir)
    copy_files_to_dir(de_paths, [de_dir for i in de_paths])
    ###Use rename.R to copy individual counts files and rename them 
    file_list =copy_de_files(base_counts, de_dir, renameDE)
    counts_default = os.path.join(outdir, "all_counts.txt")
    #de_counts_path = os.path.join(ergatis_repository, "deseq", ergatis_pid + "_differential_expression/i1/g*/all_counts")
    #print(de_counts_path)
    ###Find better way to copy the DE normalized counts file. Temp fix.
    #de_syscmd = "cp "+ de_counts_path + " " + counts_default
    #os.system(de_syscmd)


def generate_ge_report(outdir, ergatis_repository, ergatis_pid):
    """
    Input: path to input and output directory

    Output: GE files copied
    """
    # Remove redundant references and separators in path, if applicable
    ##Only run GE if HTSEQ was run!!!
    outdir = os.path.normpath(outdir)
    ergatis_repository = os.path.normpath(ergatis_repository)
    print("In GE")
    ge_paths = [os.path.join(ergatis_repository, "rpkm_coverage_stats", ergatis_pid + "_rpkm_cvg/i1/g*/genic_coverage/*.txt"), os.path.join(ergatis_repository, "htseq", ergatis_pid + "*_counts/i1/g*/*.counts")]
    ge_dir = ["RPKM","Counts"]
    ge_full_paths = prepend(ge_dir, outdir)
    print(ge_paths[1])
    makedir(ge_full_paths)
    copy_files_to_dir(ge_paths, ge_full_paths)

    # SAdkins Note - I would rather not include a RPKM or Counts output directory if that component was not run, but for now
    # this component will create the output directories and report_generator.py will check if empty or not.

    if glob.glob(ge_paths[1]):
        print("Generating all counts from htseq files")
        all_counts = generate_all_counts( ge_full_paths[1])
        counts_out = os.path.join(outdir, "all_counts.txt")
        all_counts.to_csv(path_or_buf = counts_out, header = True, sep = "\t", index = False)
    else:
        print("Unable to locate counts files.")

def generate_fastqc_report(outdir, ergatis_repository, ergatis_pid):
    """
    Input: path to input and output directory

    Output: Path and Name to FastQC report
    """
    # Remove redundant references and separators in path, if applicable
    outdir = os.path.normpath(outdir)
    ergatis_repository = os.path.normpath(ergatis_repository)
    ###This command is failing so I am changing it back to not use os.join for now
    #fqc_path = os.path.join(ergatis_repository, "/fastqc_stats/", ergatis_pid, "_fastqc/i1/g*/*/Images/")
    fqc_path = os.path.normpath(ergatis_repository+"/fastqc_stats/"+ergatis_pid+ "_fastqc/i1/g*/*/Images/")
    print(ergatis_repository) 
    print("Staring Directory Setup")
    print(fqc_path)
    dir_list = ["FastQC_Files", "FastQC_Files/KmerProfiles", "FastQC_Files/BaseQuality", "FastQC_Files/adapter_content", "FastQC_Files/duplication_levels", "FastQC_Files/per_base_n_content", "FastQC_Files/per_base_sequence_content", "FastQC_Files/per_sequence_gc_content", "FastQC_Files/per_sequence_quality", "FastQC_Files/per_tile_quality", "FastQC_Files/sequence_length_distribution"]

    extension_list = ["*_sequence.kmer_profiles.png", "*_base_quality.png", "*_sequence.adapter_content.png", "*_sequence.duplication_levels.png", "*_sequence.per_base_n_content.png", "*_sequence.per_base_sequence_content.png", "*_sequence.per_sequence_gc_content.png", "*_sequence.per_sequence_quality.png", "*_sequence.per_tile_quality.png", "*_sequence.sequence_length_distribution.png"]

    copy_to = ["FastQC_Files/KmerProfiles", "FastQC_Files/BaseQuality", "FastQC_Files/adapter_content", "FastQC_Files/duplication_levels", "FastQC_Files/per_base_n_content", "FastQC_Files/per_base_sequence_content", "FastQC_Files/per_sequence_gc_content", "FastQC_Files/per_sequence_quality", "FastQC_Files/per_tile_quality", "FastQC_Files/sequence_length_distribution"]

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
        if not os.path.isdir(path):
            os.makedirs(path)
    return True


def prepend(files, dirpath):
    """
    Input: List of files in the input directory and path to input directory	    Input: List of files in the input directory and path to input directory

    Output: List of full paths to log files in input directory	    Output: List of full paths to log files in input directory
    """
    return [os.path.join(dirpath, i) for i in files]

if __name__ == '__main__':
    main()
