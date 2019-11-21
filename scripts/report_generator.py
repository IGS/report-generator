import sys, os
import argparse
import subprocess
import pandas
from bdbag import bdbag_api
from shutil import copy

def main():
    parser = argparse.ArgumentParser( description='Report Generation Utility')
    parser.add_argument("-b", "--bdbag", dest="bdbag",help="Path to bdbag archive", metavar="BDBAG")
    parser.add_argument("-o", "--outdir", dest="outdir",help="Path to Output directory", metavar="PATH")
    parser.add_argument("-n", "--name", dest="pname",help="Name of Project must match that when generating bag", metavar="NAME")
    parser.add_argument("-1", "--fastqc", dest="fqc",help="Make fastQC Report", action='store_true')
    parser.add_argument("-2", "--align", dest="aln",help="Make Alignment Report", action='store_true')
    parser.add_argument("-3", "--ge", dest="ge",help="Make GE Report", action='store_true')
    parser.add_argument("-4", "--de", dest="de",help="Make DE Report", action='store_true')
    parser.add_argument("-a", "--all", dest="all",help="Make All Reports", action='store_true')
    parser.add_argument("-p", "--prok", dest="prok",help="Make Prok Reports", action='store_true')
    parser.add_argument("-u", "--update", dest="update",help="Update BdBag", action='store_true')
    parser.add_argument("-i", "--info", dest="info",help="Path to Info file", metavar="PATH")
    parser.add_argument("-m", "--mapping", dest="mapping",help="Path to mapping file", metavar="PATH")

    parser.add_argument("-f", "--fqcwrapper", dest="fwrap",help="Path to fastQC wrapper script", metavar="PATH")
    parser.add_argument("-w", "--alignwrapper", dest="awrap",help="Path to Alignment wrapper script", metavar="PATH")
    parser.add_argument("-d", "--dewrapper", dest="dwrap",help="Path to DE wrapper script", metavar="PATH")
    parser.add_argument("-g", "--gewrapper", dest="gwrap",help="Path to GE wrapper script", metavar="PATH")

    args = parser.parse_args()
    #wrap_FQC="/usr/local/packages/report_generation/wrapper_FastQC.R"
    if args.fwrap:
        wrap_FQC=args.fwrap
    else:
        wrap_FQC="/home/apaala.chatterjee/RNA_Report/wrapper_FastQC.R"

    if args.awrap:
        wrap_ALN = args.awrap
    elif not args.awrap and args.prok:
        wrap_ALN = "/home/apaala.chatterjee/RNA_Report/wrapper_Alignment_prok.R"
    else:
        wrap_ALN = "/home/apaala.chatterjee/RNA_Report/wrapper_Alignment.R"

    if args.gwrap:
        wrap_GE = args.gwrap
    elif not args.gwrap and args.mapping:
        wrap_GE_mapping="/home/apaala.chatterjee/RNA_Report/wrapper_GE_mapping.R"
    else:
        wrap_GE="/home/apaala.chatterjee/RNA_Report/wrapper_GE.R"

    if args.dwrap:
        wrap_DE = args.dwrap
    elif not args.dwrap and args.mapping:
        wrap_DE_mapping="/home/apaala.chatterjee/RNA_Report/wrapper_DE_mapping.R"
    else:
        wrap_DE="/home/apaala.chatterjee/RNA_Report/wrapper_DE.R"

    #wrap_ALN_prok="/home/apaala.chatterjee/RNA_Report/wrapper_Alignment_prok.R"
    #wrap_ALN= "/home/apaala.chatterjee/RNA_Report/wrapper_Alignment.R"
    #wrap_GE="/home/apaala.chatterjee/RNA_Report/wrapper_GE.R"
    #wrap_DE="/home/apaala.chatterjee/RNA_Report/wrapper_DE.R"
    #wrap_DE_mapping="/home/apaala.chatterjee/RNA_Report/wrapper_DE_mapping.R"
    #wrap_GE_mapping="/home/apaala.chatterjee/RNA_Report/wrapper_GE_mapping.R"
    counts_script="/usr/local/packages/report_generation/Generate_all_counts.R"
    rpath="/usr/local/packages/r-3.4.0/bin/Rscript"

    if args.bdbag and args.outdir and args.pname:
        extracted_path = extract_bag(args.bdbag, output_directory=args.outdir, project_name=args.pname)
        copy(args.info, extracted_path)

        if args.all:
            wrap_dir = os.path.dirname(wrap_DE)
            generate_all_reports(extracted_path, wrap_dir, rpath, args.pname, args.info, args.prok, args.mapping, args)
        else:
            if args.fqc:
                generate_fastqc_report(extracted_path, wrap_FQC, rpath, args.pname, args.info)

            if args.aln:
                generate_alignment_report(extracted_path, wrap_ALN, rpath, args.pname, args.info)

            if args.ge:
                if args.mapping:
                    print("In mapping file GE")
                    map_to_GE(extracted_path, wrap_GE_mapping, rpath, args.pname, args.info, args.mapping)
                else:
                    generate_ge_report(extracted_path, wrap_GE, rpath, args.pname, args.info)

            if args.de:
                if args.mapping:
                    print("In mapping file DE")
                    map_to_DE(extracted_path, wrap_DE_mapping, rpath, args.pname, args.info, args.mapping)
                else:
                    generate_de_report(extracted_path, wrap_DE, rpath, args.pname, args.info)

        if(args.update):
            project_outdir=os.path.join(args.outdir,args.pname)
            update_bag(project_outdir)
            #bdbag_api.make_bag(project_outdir, update=True)
            #bdbag_api.archive_bag(project_outdir, "zip")

def extract_bag(bdbag_zip_path, output_directory=None, project_name=None):
    """Extract BDBag contents into named output directory in original BDBag location."""
    (before, sep, after) = bdbag_zip_path.rpartition('.zip')
    prefix = os.path.basename(before)
    if project_name:
        prefix = project_name
    outdir = os.path.dirname(before)
    if output_directory:
        outdir = output_directory
    outdir = os.path.normpath(outdir)
    bdbag_api.extract_bag(bdbag_zip_path, output_path=outdir)
    return os.path.join(outdir, prefix, "data")

def generate_all_counts(path_to_counts):
    """
    Input: path to directory with counts files

    Output: Single merged counts dataframe
    """

    counts_list = [f for f in os.listdir(path_to_counts) if f.endswith('.counts') or f.endswith('.counts.txt')]
    appended_counts_list = [os.path.join(path_to_counts, counts_list) for i in counts_list]
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

def generate_all_reports(extracted_path, wrappers_dir, rpath, project_name, info_file, prok, mapping_file, args=None):
    """Generate all possible reports.  Assumes all wrapper scripts are in the same directory."""

    wrap_FQC = os.path.join(wrappers_dir, "wrapper_FastQC.R")
    if args.fwrap:
        wrap_FQC = args.fwrap
    generate_fastqc_report(extracted_path, wrap_FQC, rpath, args.pname, args.info)

    wrap_ALN = os.path.join(wrappers_dir, "wrapper_Alignment.R")
    if prok:
        wrap_ALN_prok = os.path.join(wrappers_dir, "wrapper_Alignment_prok.R")
    if args.awrap:
        wrap_ALN = args.awrap
    generate_alignment_report(extracted_path, wrap_ALN, rpath, args.pname, args.info)

    if mapping_file:
        wrap_GE_mapping = os.path.join(wrappers_dir, "wrapper_GE_mapping.R")
        wrap_DE_mapping = os.path.join(wrappers_dir, "wrapper_DE_mapping.R")
        if args.gwrap:
            wrap_GE_mapping = args.gwrap
        if args.dwrap:
            wrap_DE_mapping = args.dwrap
        map_to_GE(extracted_path, wrap_GE_mapping, rpath, args.pname, args.info, args.mapping)
        map_to_DE(extracted_path, wrap_DE_mapping, rpath, args.pname, args.info, args.mapping)
    else:
        wrap_GE = os.path.join(wrappers_dir, "wrapper_GE.R")
        wrap_DE = os.path.join(wrappers_dir, "wrapper_DE.R")
        if args.gwrap:
            wrap_GE = args.gwrap
        if args.dwrap:
            wrap_DE = args.dwrap
        generate_ge_report(extracted_path, wrap_GE, rpath, args.pname, args.info)
        generate_de_report(extracted_path, wrap_DE, rpath, args.pname, args.info)

def generate_alignment_report(outdir, wrapper, rpath, project_name, info_file):
    """
    Input: path to output directory with extracted files, path to wrapper R script, path to R, name of project and path to info file.

    Output: Alignment report generated and output plots saved in AlignmentFiles
    """
    outdir = os.path.normpath(outdir)
    output_files = os.path.join(outdir, "AlignmentFiles")
    if not os.path.isdir(output_files):
        os.makedirs(output_files)
    try:
        subprocess.check_call([rpath, wrapper, project_name, outdir, info_file, outdir], shell = False)
    except subprocess.CalledProcessError as e:
        print(e)

def generate_de_report(outdir, wrapper, rpath, project_name, info_file):
    """
    Input: path to output directory with extracted files, path to wrapper R script, path to R, name of project and path to info file.

    Output: DE report generated and outputs saved in DE_Outputs
    """

    outdir = os.path.normpath(outdir)
    output_files = os.path.join(outdir, "DE_Outputs")
    all_count_path = os.path.join(outdir, "all_counts.txt")
    if not os.path.isdir(output_files):
        os.makedirs(output_files)
    try:
        subprocess.check_call([rpath, wrapper, project_name, outdir, info_file, outdir], shell = False)
    except subprocess.CalledProcessError as e:
        print(e)

def generate_fastqc_report(outdir, wrapper, rpath, project_name, info_file):
    """
    Input: path to output directory with extracted files, path to wrapper R script, path to R, name of project and path to info file.

    Output: FastQC report generate and output files FastQC_Outputs
    """
    outdir = os.path.normpath(outdir)
    output_files = os.path.join(outdir, "FastQC_Outputs")
    if not os.path.isdir(output_files):
        os.makedirs(output_files)
    try:
        subprocess.check_call([rpath, wrapper, project_name, outdir, info_file, outdir], shell = False)
    except subprocess.CalledProcessError as e:
        print(e)
    #syscmd=rpath +" "+ wrap_FQC +" "+ args.pname +" "+ pdir +" " + args.info + " "+pdir

def generate_ge_report(outdir, wrapper, rpath, project_name, info_file):
    """
    Input: path to output directory with extracted files, path to wrapper R script, path to R, name of project and path to info file.

    Output: GE report generate and output files in GE_Outputs
    """

    outdir = os.path.normpath(outdir)
    output_files = os.path.join(outdir, "GE_Outputs")
    all_count_path = os.path.join(outdir, "all_counts.txt")
    count_file_expected_path = os.path.join(outdir, "Counts")
    if not os.path.isdir(output_files):
        os.makedirs(output_files)

    # Ensure counts directory has files in BDBag (AKA the 'htseq' component was in the pipeline)
    if os.listdir(count_file_expected_path):
        if not os.path.exists(all_count_path):
            merged_counts = generate_all_counts(count_file_expected_path)
            merged_counts.to_csv(path_or_buf = all_count_path, header = True, sep = "\t", index = False)
        try:
            subprocess.check_call([rpath, wrapper, project_name, outdir, info_file, outdir], shell = False)
        except subprocess.CalledProcessError as e:
            print(e)

def map_to_DE(outdir, wrapper, rpath, project_name, info_file, mappingf):
    """
    Input: path to output directory with extracted files, path to wrapper R script, path to R, name of project and path to info file.

    Output: DE report generate and output files in DE_Outputs
    """
    outdir = os.path.normpath(outdir)
    output_files = os.path.join(outdir, "DE_Outputs")
    all_count_path = os.path.join(outdir, "all_counts.txt")
    if not os.path.isdir(output_files):
        os.makedirs(output_files)
    try:
        subprocess.check_call([rpath, wrapper, project_name, outdir, info_file, outdir, mappingf], shell = False)
    except subprocess.CalledProcessError as e:
        print(e)

def map_to_GE(outdir, wrapper, rpath, project_name, info_file, mappingf):
    """
    Input: path to output directory with extracted files, path to wrapper R script, path to R, name of project, path to info file and path to mapping file

    Output: GE report generate and output files in GE_Outputs
    """

    outdir = os.path.normpath(outdir)
    output_files = os.path.join(outdir, "GE_Outputs")
    all_count_path = os.path.join(outdir, "all_counts.txt")
    if not os.path.isdir(output_files):
        os.makedirs(output_files)
    try:
        subprocess.check_call([rpath, wrapper, project_name, outdir, info_file, outdir, mappingf], shell = False)
    except subprocess.CalledProcessError as e:
        print(e)

def update_bag(outdir):
    bdbag_api.make_bag(outdir, update=True)
    return bdbag_api.archive_bag(outdir, "zip")

if __name__ == '__main__':
    main()
