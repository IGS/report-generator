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

    (options, args) = parser.parse_args()
    #wrap_FQC="/usr/local/packages/report_generation/wrapper_FastQC.R"
    if options.fwrap:
        wrap_FQC=options.fwrap
    else:
        wrap_FQC="/home/apaala.chatterjee/RNA_Report/wrapper_FastQC.R"
    if options.awrap:
        wrap_ALN = options.awrap
    elif not options.awrap and options.prok:
        wrap_ALN_prok="/home/apaala.chatterjee/RNA_Report/wrapper_Alignment_prok.R"
    else:
        wrap_ALN= "/home/apaala.chatterjee/RNA_Report/wrapper_Alignment.R"
    if options.gwrap:
        wrap_GE = options.gwrap
    elif not options.gwrap and options.mapping:
        wrap_GE_mapping="/home/apaala.chatterjee/RNA_Report/wrapper_GE_mapping.R"
    else:
        wrap_GE="/home/apaala.chatterjee/RNA_Report/wrapper_GE.R"
    if options.dwrap:
        wrap_DE = options.dwrap
    elif not options.dwrap and options.mapping:
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

    if options.bdbag and options.outdir and options.pname:
        extracted_path = extract_bag(options.bdbag, output_directory=options.outdir, project_name=options.pname)
        copy(options.info, extracted_path)

        if options.all:
            wrap_dir = os.path.dirname(wrap_DE)
            generate_all_reports(extracted_path, wrap_dir, rpath, options.pname, options.info, options.prok, options.mapping, options)
        else:
            if options.fqc:
                generate_fastqc_report(extracted_path, wrap_FQC, rpath, options.pname, options.info)

            if options.aln:
                if options.prok:
                    generate_alignment_report_prok(extracted_path, wrap_ALN_prok, rpath, options.pname, options.info)
                else:
                    generate_alignment_report(extracted_path, wrap_ALN, rpath, options.pname, options.info)

            if options.ge:
                if options.mapping:
                    print("In mapping file GE")
                    map_to_GE(extracted_path, wrap_GE_mapping, rpath, options.pname, options.info, options.mapping)
                else:
                    generate_ge_report(extracted_path, wrap_GE, rpath, options.pname, options.info)

            if options.de:
                if options.mapping:
                    print("In mapping file DE")
                    map_to_DE(extracted_path, wrap_DE_mapping, rpath, options.pname, options.info, options.mapping)
                else:
                    generate_de_report(extracted_path, wrap_DE, rpath, options.pname, options.info)

        if(options.update):
            update_bag(options.outdir, project_name=options.pname)

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
    return os.path.normpath(outdir + "/" + prefix + "/data/")

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

def generate_all_reports(outdir, wrappers_dir, rpath, project_name, info_file, prok, mapping_file, opts=None):
    """Generate all possible reports.  Assumes all wrapper scripts are in the same directory."""

    wrap_FQC = os.path.join(wrappers_dir, "wrapper_FastQC.R")
    if opts.fwrap:
        wrap_FQC = opts.fwrap
    generate_fastqc_report(extracted_path, wrap_FQC, rpath, options.pname, options.info)

    if prok:
        wrap_ALN_prok = os.path.join(wrappers_dir, "wrapper_Alignment_prok.R")
        if opts.awrap:
            wrap_ALN_prok = opts.awrap
        generate_alignment_report_prok(extracted_path, wrap_ALN_prok, rpath, options.pname, options.info)
    else:
        wrap_ALN = os.path.join(wrappers_dir, "wrapper_Alignment.R")
        if opts.awrap:
            wrap_ALN = opts.awrap
        generate_alignment_report(extracted_path, wrap_ALN, rpath, options.pname, options.info)

    if mapping_file:
        wrap_GE_mapping = os.path.join(wrappers_dir, "wrapper_GE_mapping.R")
        wrap_DE_mapping = os.path.join(wrappers_dir, "wrapper_DE_mapping.R")
        if opts.gwrap:
            wrap_GE_mapping = opts.gwrap
        if opts.dwrap:
            wrap_DE_mapping = opts.dwrap
        map_to_GE(extracted_path, wrap_GE_mapping, rpath, options.pname, options.info, options.mapping)
        map_to_DE(extracted_path, wrap_DE_mapping, rpath, options.pname, options.info, options.mapping)
    else:
        wrap_GE = os.path.join(wrappers_dir, "wrapper_GE.R")
        wrap_DE = os.path.join(wrappers_dir, "wrapper_DE.R")
        if opts.gwrap:
            wrap_GE = opts.gwrap
        if opts.dwrap:
            wrap_DE = opts.dwrap
        generate_ge_report(extracted_path, wrap_GE, rpath, options.pname, options.info)
        generate_de_report(extracted_path, wrap_DE, rpath, options.pname, options.info)

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
        subprocess.check_call([rpath, wrapper, project_name, outdir, info_file, outdir], shell = False)
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
        subprocess.check_call([rpath, wrapper, project_name, outdir, info_file, outdir], shell = False)
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

    output_files = os.path.normpath(outdir + "/DE_Outputs/")
    if not os.path.exists(output_files):
        os.makedirs(output_files)
    all_count_path = os.path.normpath(outdir + "/all_counts.txt")
    try:
        subprocess.check_call([rpath, wrapper, project_name, outdir, info_file, outdir, mappingf], shell = False)
    except subprocess.CalledProcessError as e:
        print(e)

def map_to_GE(outdir, wrapper, rpath, project_name, info_file, mappingf):
    """
    Input: path to output directory with extracted files, path to wrapper R script, path to R, name of project, path to info file and path to mapping file

    Output: GE report generate and output files in GE_Outputs
    """

    output_files = os.path.normpath(outdir + "/GE_Outputs/")
    if not os.path.exists(output_files):
        os.makedirs(output_files)
    all_count_path = os.path.normpath(outdir + "/all_counts.txt")
    try:
        subprocess.check_call([rpath, wrapper, project_name, outdir, info_file, outdir, mappingf], shell = False)
    except subprocess.CalledProcessError as e:
        print(e)

def update_bag(outdir):
    bdbag_api.make_bag(outdir, update=True)
    return bdbag_api.archive_bag(outdir, "zip")

if __name__ == '__main__':
    main()
