<%@ page language="java" contentType="text/html; charset=UTF-8" pageEncoding="UTF-8"%>
<%@ taglib prefix="s" uri="/struts-tags"%>
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
    <!-- Meta, title, CSS, favicons, etc. -->
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <link rel="shortcut icon" href="../images/favicon.png"> 
    <title>CFEA</title>

    <!-- Bootstrap -->
    <link href="../css/bootstrap.min.css" rel="stylesheet">
    <!-- Font Awesome -->
    <link href="../css/font-awesome/css/font-awesome.min.css" rel="stylesheet">
    <!-- NProgress -->
    <link href="../css/nprogress.css" rel="stylesheet">
    <!-- Custom Theme Style -->
    <link href="../css/custom.min.css" rel="stylesheet">
    <link rel="stylesheet" href="../css/loading.css" type="text/css" media="screen" />
		<!-- prettify print -->
		<script type="text/javascript" src="../js/prettify.js"></script>
		<link rel="stylesheet" type="text/css" href="../css/prettify.css"/>
    <style>
    	.help{
    		font-size:16px;
    	}
    	.content{
    		font-size:16px;
    		margin-bottom:20px;
    		border: 1px solid rgb(238,232,225);
    	}
    </style>
  </head>

  <body class="nav-md" onload="prettyPrint()">
    <div class="container body">
      <div class="main_container">
        <div class="col-md-3 left_col">
          <div class="left_col scroll-view">
            <div class="navbar nav_title" style="border: 0;">
              <a href="/CFEA/index.jsp" class="site_title"><span>CFEA</span></a>
            </div>
            <div class="clearfix"></div>
            <!-- /menu profile 450K info -->
            <br />
            <!-- sidebar menu -->
            <div id="sidebar-menu" class="main_menu_side hidden-print main_menu">
              <div class="menu_section">
                <h3>Welcome to CFEA</h3>
                <ul class="nav side-menu">
                  <li><a href="/CFEA/index.jsp"><i class="fa fa-home"></i> Home <span class="fa fa-chevron-right"></span></a></li>
                  <li><a href="/CFEA/jsp/Cancer-stage.jsp"><i class="fa fa-film"></i> Cancer-stage browser <span class="fa fa-chevron-right"></span></a></li>
                  <li><a href="/CFEA/jsp/Pipelines.jsp"><i class="fa fa-film"></i> CFEA pipelines <span class="fa fa-chevron-right"></span></a></li>
                  <li><a href="/CFEA/jsp/Statistics.jsp"><i class="fa fa-film"></i> Statistics <span class="fa fa-chevron-right"></span></a></li>
                  <li><a href="/CFEA/jsp/Tutorial.jsp"><i class="fa fa-book"></i> Tutorial <span class="fa fa-chevron-right"></span></a></li>
                  <li><a href="/CFEA/jsp/Download.jsp"><i class="fa fa-download"></i> Download <span class="fa fa-chevron-right"></span></a></li>
                 <!--  <li><a href="/MEddc/jsp/Contact.jsp"><i class="fa fa-envelope"></i> Contact <span class="fa fa-chevron-right"></span></a></li></li> -->
                </ul>
              </div>
            </div>
          </div>
        </div>
        <!-- top navigation -->
<!--         <div class="top_nav"> -->
<!--           <div class="nav_menu"> -->
<!--             <nav style="height:60px;"> -->
<!--               <div class="nav toggle"> -->
<!--                 <a id="menu_toggle"><i class="fa fa-bars"></i></a> -->
<!--               </div> -->
<!--             </nav> -->
<!--           </div> -->
<!--         </div> -->
        <!-- /top navigation -->
        <!-- page content -->
        <div class="fuzzy-advSearch" id="fuzzy"></div>
	      <div class="loading" id="loading">
		  <img alt="loading" title="loading" src="../images/index.gif">
        </div>
        <div class="right_col" role="main">
          <div class="">
            <div class="row">
              <div class="col-md-12 col-sm-12 col-xs-12">
                  <div class="x_panel">
                  <div class="x_title">
                  	<h2><i class="fa fa-book"></i> Pipelines</h2>
                    <ul class="nav navbar-right panel_toolbox">
                      <li><a class="collapse-link"><i class="fa fa-chevron-up"></i></a>
                      </li>
                    </ul>
                    <div class="clearfix"></div>
                  </div>
                  <div class="x_content">
<!--                   <br/> -->
                    <div class="content">
                       <div>
                       	<B style="margin:5px 5px 5px 5px;">CFEA pipelines</B>
                        <form style="margin:15px 15px 15px 35px;">
                        	<ul>
	                        	<li><a href="#Overview">Overview of Pipelines</a></li>
	                        	<li><a href="#5HMC">5HMC</a></li>
                            <ul>
	                        			<li><a href="#callpeak">callpeak.py</a></li>
	                        			<li><a href="#mapping_bio_13">mapping_bio_13.py</a></li>
	                        		</ul>
	                        	<li><a href="#450K">450K</a></li>
														<ul>
	                        			<li><a href="#parse_450K">parse_450K.R</a></li>
	                        		</ul>
	                        	<li><a href="#general_step">General_step</a></li>
                            <ul>
	                        			<li><a href="#beta2wig">beta2wig.py</a></li>
	                        			<li><a href="#cov2wig">cov2wig.py</a></li>
                                <li><a href="#fastqc0">fastqc0.py</a></li>
	                        			<li><a href="#get_reads_samples_multi">get_reads_samples_multi.py</a></li>
	                        		</ul>
	                        	<li><a href="#MCTA-seq">MCTA-seq</a>
	                        		<ul>
	                        			<li><a href="#ch3deleate">ch3deleate.pl</a></li>
	                        			<li><a href="#QCnolane">QCnolane.pl</a></li>
                                <li><a href="#rm_firstx_leny">rm_firstx_leny.pl</a></li>
	                        			<li><a href="#trim_mapping_bio_19">trim_mapping_bio_19.py</a></li>
	                        		</ul>
	                        	</li>
	                        	<li><a href="#Nucelosome">Nucelosome</a></li>
                            <ul>
	                        			<li><a href="#atacCallPeak">atacCallPeak.sh</a></li>
	                        			<li><a href="#atacMapTrim">atacMapTrim.sh</a></li>
	                        		</ul>
	                        	<li><a href="#RRBS_WGBS">RRBS_WGBS</a></li>
                            <ul>
	                        			<li><a href="#extract_cov_li_20">extract_cov_li_20.py</a></li>
                                <li><a href="#get_meta_sql_li_20">get_meta_sql_li_20.py</a></li>
                                <li><a href="#getCG-1.0.0">getCG-1.0.0.pl</a></li>
	                        			<li><a href="#grep_cpg_li_20">grep_cpg_li_20.py</a></li>
                                <li><a href="#mergeFq">mergeFq.py</a></li>
	                        			<li><a href="#trim_mapping_li_20">trim_mapping_li_20.py</a></li>
	                        		</ul>
                        	</ul>
                        </form> 
                       </div>              	
                    </div>
                    <div class="help">
	                     <div id="Overview" style="margin-top:10px;border:1px solid rgb(238,232,225);">
	                     	<B style="margin:5px 5px 5px 5px;">Overview of Pipelines</B><br/>
		                     <p style="margin:15px 15px 15px 15px;">
		                     	MutExGenome is a resource that comprehensive curated mutually exclusive alterations in cancer genome and integrated Fisher’s statistically inferred from The Cancer Genome Atlas. Currently, MutExGenome documents 75,714 mutually exclusive relationships involving 140 human cancers of 27 tissue origins. MutExGenome provides a user-friendly interface for conveniently browsing pan-cancer map of oncogenic dependencies, and mutually exclusive network of genetic alterations, searching and downloading the mutually exclusive genomic events in various cancers. In addition, MutExGenome provides a submission function that allows researchers to submit novel mutually exclusive associations that are not documented. 
								<br/>The detailed usage of the database is as followings:
		                     </p>
	                     </div>
						 <div id="5HMC" style="margin-top:10px;border:1px solid rgb(238,232,225);">
						 	<B style="margin:5px 5px 5px 5px;">5HMC</B>
		                     <p style="margin:15px 15px 15px 15px;">
		                     	The home page contains two schematic diagrams to help users understand our database. The left figure is a conceptual graph illustrating the mutually exclusive genomic relationships across patients. The right figure in the home page shows a human body, some common cancers are marked as abbreviations to facilitate the user to 450Kly browse mutually exclusive pairs in the marked cancer type. For example, after clicking the 'OV' link, the search engine will run and return the results containing a table showing all related mutually exclusive events and corresponding mutually exclusive network associated with ovarian cancer(as figure below shows).
		                     </p>
                         <p style="margin:15px 15px 15px 15px;">
		                        <B ><span id="callpeak">callpeak.py</span></B><br/><!--代码位置 -->
		                     	<pre class="prettyprint linenums">
import os
import glob
import argparse
import subprocess
import multiprocessing


"""
this script is written to call peaks from the deduplicated bam files with macs2 software
__author__ = likai
__email__ = likai@wibe.ac.cn
"""


def get_bam_files(mapping_dir):
    """
    get all the deduplicated bam file from the mapping directory if there are any
    :param mapping_dir: directory to the deduplicated bam files
    :return: python list with all targeted bam files
    """
    print(f"Start reading from {mapping_dir}")
    bam_files = glob.glob(f"{mapping_dir}/*.dedup.bam")
    if not bam_files:
        bam_files = glob.glob(f"{mapping_dir}/*.sort.bam")
    return bam_files


def call_peak(bam, peak_dir, log_dir):
    """
    we use macs2 to call peak from deduplicated bam files
    :param bam: path to each deduplicated bam files
    :param peak_dir: directory to put results of macs2
    :param log_dir: directory to put log files in each step
    :return: none
    """
    prefix = os.path.basename(bam).split(".")[0]
    print(f"    Calling Peak:  {prefix}")
    this_sample_peak_out = f"{peak_dir}/{prefix}"
    peak_log = f"{log_dir}/{prefix}.macs2.log"
    error_log = f"{log_dir}/{prefix}.macs2.error.log"
    if not os.path.exists(this_sample_peak_out):
        os.makedirs(this_sample_peak_out)
    cmd_macs2 = f"macs2 callpeak -t {bam} -g hs --nomodel " \
                f"--extsize 200 -n {prefix} --outdir {this_sample_peak_out} --call-summits > {peak_log} 2>&1"
    try:
        subprocess.check_output(cmd_macs2, shell=True)
    except Exception as e:
        print(e)
        with open(error_log, "w+") as out:
            out.write(str(e) + "\n")


def multi_process_run(bam_files, num_process, peak_dir, log_dir):
    """
    to run faster, we create a processing pool with a certain number of processes, the params
    are the same to function call_peak
    :param bam_files:
    :param num_process:
    :param peak_dir:
    :param log_dir:
    :return:
    """
    pool = multiprocessing.Pool(processes=int(num_process))
    for bam in bam_files:
        pool.apply_async(call_peak, (bam, peak_dir, log_dir,))
    pool.close()
    pool.join()


def get_multiqc_files(peak_dir, peak_out_dir):
    """
    we use multiQC software to collect the qc results for each peak files by macs2
    :param peak_dir: directory to put results of macs2
    :param peak_out_dir:
    :return:
    """
    print(f"Start reading from {peak_dir}")
    peak_files = glob.glob(f"{peak_dir}/*/*_peaks.xls")
    peak_qc_path = f"{peak_out_dir}/peaks_multiqc_sample.txt"
    with open(peak_qc_path, "w+") as out:
        out.write("\n".join(peak_files) + "\n")
    return peak_qc_path


def multiqc(peak_qc_path, peak_out_dir, n):
    """
    run multiQC to collect statics for peak files of each deduplicated bam files by macs2
    :param peak_qc_path: directory to the peak file called by macs2
    :param peak_out_dir: directory to put multiqc result
    :param n:
    :return:
    """
    print("Start multiqc")
    n_peak = f"{n}_peak"
    cmd_peak_multiqc = f"multiqc -n {n_peak} -s -o {peak_out_dir} -l {peak_qc_path}"
    try:
        subprocess.check_output(cmd_peak_multiqc, shell=True)
    except Exception as e:
        print(e)


def main():
    parser = argparse.ArgumentParser(description="A program to call peaks from 5hmc bam files")
    parser.add_argument("-i", action="store", dest="mapping_dir")
    parser.add_argument("-p", action="store", dest="num_process", help="processes to use in the program", default=10)
    parser.add_argument("-m", action="store", dest="peak_dir", help="path to put mapping result in")
    parser.add_argument("-l", action="store", dest="log_dir", help="path to put log files")
    parser.add_argument("-n", action="store", dest="n", help="project name to use in MultiQC output[bioproject.13]")
    parser.add_argument("-o", action="store", dest="peak_out_dir", help="directory to put peak qc results")
    results = parser.parse_args()
    if not all([results.mapping_dir, results.num_process,
                results.peak_dir, results.log_dir, results.n, results.peak_out_dir]):
        print("too few arguments, type -h for more information")
        exit(-1)
    mapping_dir = results.mapping_dir if not results.mapping_dir.endswith("/") else results.mapping_dir.rstrip("/")
    num_process = results.num_process
    peak_dir = results.peak_dir if not results.peak_dir.endswith("/") else results.peak_dir.rstrip("/")
    log_dir = results.log_dir if not results.log_dir.endswith("/") else results.log_dir.rstrip("/")
    n = results.n
    peak_out_dir = results.peak_out_dir if not results.peak_out_dir.endswith("/") else results.peak_out_dir.rstrip("/")
    for each_dir in [peak_dir, log_dir, peak_out_dir]:
        if not os.path.exists(each_dir):
            os.makedirs(each_dir)
    bam_files = get_bam_files(mapping_dir)
    multi_process_run(bam_files, num_process, peak_dir, log_dir)
    peak_qc_path = get_multiqc_files(peak_dir, peak_out_dir)
    multiqc(peak_qc_path, peak_out_dir, n)


if __name__ == '__main__':
    main()
</pre>
		                     </p>
                         <p style="margin:15px 15px 15px 15px;">
		                        <B ><span id="mapping_bio_13">mapping_bio_13.py</span></B><br/><!--代码位置 -->
		                     	<pre class="prettyprint linenums">
import os
import glob
import argparse
import subprocess
import multiprocessing

"""
This script is developed for mapping of 5hmc data, be sure that picard
samtools, fastqc and multiqc are in your local path, one way is to use
conda. 
conda install -n py2 picard, samtools, fastqc, multiqc
source activate /share/pub/lik/soft/miniconda3/envs/py2
__author__ = likai
__email__ = likai@wibe.ac.cn
"""


def get_files(path_fastq):
    """
    get all raw fastq files
    :param path_fastq: directory to raw fastq files
    :return: list with all fastq files
    """
    print(f"Start reading from {path_fastq}")
    fastq_files = glob.glob(f"{path_fastq}/*.fastq")
    return fastq_files


def mapping_bowtie2(fastq, bowtie2_path, path_index, mapping_out_dir, bam_dir, log_dir):
    """
    we first align the fastq files with bowtie2 aligner and then remove duplicates with
    picard tool set, after that, we rerun fastqc to get the quality of the deduplicated data
    :param fastq: path to each fastq file
    :param bowtie2_path: path to the bowtie2 aligner
    :param path_index: path to index files created by bowtie2 aligner
    :param mapping_out_dir: path to put the qc results of aligned bam files
    :param bam_dir: path to put aligned bam files
    :param log_dir: path to put log files in each step
    :return: none
    """
    prefix = os.path.basename(fastq).split(".")[0]
    fastq_file_1 = f"{os.path.dirname(fastq)}/{prefix}.sra_1.fastq"
    fastq_file_2 = f"{os.path.dirname(fastq)}/{prefix}.sra_2.fastq"
    sam_file = f"{bam_dir}/{prefix}.sam"
    this_sample_log = f"{log_dir}/{prefix}.mapping.error.log"
    this_sample_picard_log = f"{log_dir}/{prefix}.picard.log"
    bowtie2_out_result = f"{mapping_out_dir}/bowtie2_out/{prefix}.txt"
    if not os.path.exists(f"{mapping_out_dir}/bowtie2_out"):
        os.makedirs(f"{mapping_out_dir}/bowtie2_out")
    cmd_bowtie2 = f"{bowtie2_path} -p 30 --end-to-end -q " \
                  f"-x {path_index} -1 {fastq_file_1} -2 {fastq_file_2} -S {sam_file} > {bowtie2_out_result} 2>&1"
    sort_bam_file = f"{bam_dir}/{prefix}.sort.bam"
    cmd_sort_samtools = f"samtools view -h {sam_file} | samtools sort -@ 30 -o {sort_bam_file}"
    dedup_bam_file = f"{bam_dir}/{prefix}.dedup.bam"
    metrics_picard_file = f"{mapping_out_dir}/picard_out/{prefix}.dedup.txt"
    if not os.path.exists(f"{mapping_out_dir}/picard_out"):
        os.makedirs(f"{mapping_out_dir}/picard_out")
    cmd_picard_dedup = f"picard MarkDuplicates INPUT={sort_bam_file} OUTPUT={dedup_bam_file} " \
                       f"METRICS_FILE={metrics_picard_file} VALIDATION_STRINGENCY=LENIENT " \
                       f"ASSUME_SORTED=true REMOVE_DUPLICATES=true > {this_sample_picard_log} 2>&1"
    fastqc_dedup_out_dir = f"{mapping_out_dir}/fastqc_dedup/{prefix}"
    if not os.path.exists(fastqc_dedup_out_dir):
        os.makedirs(fastqc_dedup_out_dir)
    cmd_fastqc_dedup_bam = f"fastqc -q --extract -o {fastqc_dedup_out_dir} {dedup_bam_file}" # run fastqc for bam files
    try:
        print(f"  Aligning: {prefix}")
        subprocess.check_output(cmd_bowtie2, shell=True)
        print(f"  Sorting: {prefix}")
        subprocess.check_output(cmd_sort_samtools, shell=True)
        print(f"  Removing duplicates: {prefix}")
        subprocess.check_output(cmd_picard_dedup, shell=True)
        print(f"  Dedup QC: {prefix}")
        subprocess.check_output(cmd_fastqc_dedup_bam, shell=True)
    except Exception as e:
        print(e)
        with open(this_sample_log, "w+") as out:
            out.write(str(e) + "\n")
    

def multi_process_run(fastq_files, num_process, bowtie2_path, path_index, mapping_out_dir, bam_dir, log_dir):
    """
    to run faster, we create a processing pool with a certain number of processes, the params
    are the same to function mapping_bowtie2
    :param fastq_files:
    :param num_process:
    :param bowtie2_path:
    :param path_index:
    :param mapping_out_dir:
    :param bam_dir:
    :param log_dir:
    :return: none
    """
    pool = multiprocessing.Pool(processes=int(num_process))
    for fastq in fastq_files:
        if not fastq.endswith(".sra_2.fastq"):
            pool.apply_async(mapping_bowtie2, (fastq, bowtie2_path, path_index, mapping_out_dir, bam_dir, log_dir,))
    pool.close()
    pool.join()


def get_multiqc_file(mapping_out_dir):
    """
    we use multiQC software to collect the qc results for the aligned bam files and deduplicated bam files
    this function is designed to get all the aligned and deduplicated bam files from the mapping directory
    :param mapping_out_dir: directory to the qc results of the aligned bam files and deduplicated bam files
    :return: path to multiQC result file of aligned bam files and deduplicated bam files
    """
    print(f"Start reading from {mapping_out_dir}")
    bowtie2_out_result_dir = f"{mapping_out_dir}/bowtie2_out"
    bowtie2_out_files = glob.glob(f"{bowtie2_out_result_dir}/*.txt")
    bowtie2_out_multiqc_path = f"{mapping_out_dir}/bowtie2_multiqc_sample.txt"
    dedup_qc_dir = f"{mapping_out_dir}/fastqc_dedup"
    dedup_qc_zip_files = glob.glob(f"{dedup_qc_dir}/*/*.zip")
    dedup_qc_multiqc_path = f"{mapping_out_dir}/dedup_multiqc_sample.txt"
    with open(bowtie2_out_multiqc_path, "w+") as out:
        out.write("\n".join(bowtie2_out_files) + "\n")
    with open(dedup_qc_multiqc_path, "w+") as out:
        out.write("\n".join(dedup_qc_zip_files) + "\n")
    return bowtie2_out_multiqc_path, dedup_qc_multiqc_path


def multiqc(bowtie2_out_multiqc_path, dedup_qc_multiqc_path, n, mapping_out_dir):
    """
    run multiQC for the deduplicated bam files by picard
    :param bowtie2_out_multiqc_path:
    :param dedup_qc_multiqc_path:
    :param n:
    :param mapping_out_dir:
    :return:
    """
    print("Start multiqc")
    n_bowtie2 = f"{n}_bowtie2"
    # run multiQC for aligned bam files
    cmd_bowtie2_multiqc = f"multiqc -n {n_bowtie2} -s -o {mapping_out_dir} -l {bowtie2_out_multiqc_path}"
    n_dedup_qc = f"{n}_dedup_qc"
    # run multiQC for deduplicated bam files
    cmd_dedup_multiqc = f"multiqc -n {n_dedup_qc} -s -o {mapping_out_dir} -l {dedup_qc_multiqc_path}"
    try:
        subprocess.check_output(cmd_bowtie2_multiqc, shell=True)
        subprocess.check_output(cmd_dedup_multiqc, shell=True)
    except Exception as e:
        print(e)


def main():
    parser = argparse.ArgumentParser(description="A program to map raw fastq and remove duplicates")
    parser.add_argument("-i", action="store", dest="path_fastq")
    parser.add_argument("-p", action="store", dest="num_process", help="processes to use in the program", default=10)
    parser.add_argument("-m", action="store", dest="mapping_out_dir", help="path to put mapping result in")
    parser.add_argument("-l", action="store", dest="log_dir", help="path to put log files")
    parser.add_argument("-b", action="store", dest="bam_dir", help="path to put bam files")
    parser.add_argument("-n", action="store", dest="n", help="project name to use in MultiQC output[bioproject.13]")
    parser.add_argument("-t", action="store", dest="bowtie2_path", help="path to bowtie2 executable")
    parser.add_argument("-d", action="store", dest="path_index", help="path to bowtie2 indexes")
    results = parser.parse_args()
    if not all([results.path_fastq, results.num_process, results.mapping_out_dir, results.log_dir, 
                results.n, results.bam_dir, results.bowtie2_path, results.path_index]):
        print("too few arguments, type -h for more information!")
        exit(-1)
    path_fastq = results.path_fastq if not results.path_fastq.endswith("/") else results.path_fastq.rstrip("/")
    num_process = results.num_process
    mapping_out_dir = results.mapping_out_dir if not results.mapping_out_dir.endswith("/") else results.mapping_out_dir.rstrip("/")
    log_dir = results.log_dir if not results.log_dir.endswith("/") else results.log_dir.rstrip("/")
    bam_dir = results.bam_dir if not results.bam_dir.endswith("/") else results.bam_dir.rstrip("/")
    n = results.n
    bowtie2_path = results.bowtie2_path if not results.bowtie2_path.endswith("/") else results.bowtie2_path.rstrip("/")
    path_index = results.path_index if not results.path_index.endswith("/") else results.path_index.rstrip("/")
    for each_dir in [mapping_out_dir, log_dir, bam_dir]:
        if not os.path.exists(each_dir):
            os.makedirs(each_dir)
    fastq_files = get_files(path_fastq)
    multi_process_run(fastq_files, num_process, bowtie2_path, path_index, mapping_out_dir, bam_dir, log_dir)
    bowtie2_out_multiqc_path, dedup_qc_multiqc_path = get_multiqc_file(mapping_out_dir)
    multiqc(bowtie2_out_multiqc_path, dedup_qc_multiqc_path, n, mapping_out_dir)

    
if __name__ == '__main__':
    main()
</pre>
		                     </p>
						 </div>
	                     
	                     <div id="450K" style="margin-top:10px;border:1px solid rgb(238,232,225);">
	                     	<B style="margin:5px 5px 5px 5px;">450K</B>
		                     <p style="margin:15px 15px 15px 15px;"><!--代码位置 -->
		                     	On the home page, users also can search mutually exclusive events by inputting the gene symbol or cancer name of interest. In this page, only one term could be used, for combined search, please go to “Search” module. For example, after inputting BRCA1, the search engine will run and return the results containing a table showing all mutually exclusive partners associated with BRCA1 and corresponding mutually exclusive network (as figure below shows).
		                     </p>
												 <p style="margin:15px 15px 15px 15px;">
		                        <B ><span id="parse_450K">parse_450K.R</span></B><br/><!--代码位置 -->
		                     	<pre class="prettyprint linenums">
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(missMethyl)
library(matrixStats)
library(minfiData)
library(readr)


myArgs <- commandArgs(trailingOnly = T)
# path to your idat files
path <- myArgs[1]
plot1 <- myArgs[2]
plot2 <- myArgs[3]
pathBValue <- myArgs[4]

if(length(myArgs) != 4){
    stop("Usage: Rscript [path_to_450K_idat] [p_values.pdf] [density_plot.pdf] [path_to_beta_value]")
}

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")

cat("Start reading in ", path, "\n")
rgSet <- read.metharray.exp(base = path, recursive = F, force=TRUE)
# rgSet <- read.metharray.exp(base = path, recursive = F)
sampleName <- sapply(colnames(rgSet), FUN = function(X){unlist(strsplit(X, split = "_"))[1]})
colnames(rgSet) <- sampleName

# detect P value
detP <- detectionP(rgSet)
pdf(plot1)
barplot(colMeans(detP), col = "steelblue", las = 2, cex.names = 0.8, ylab = "Mean detection p-values")
dev.off()

keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]
detP <- detP[,keep]
mSetSq <- preprocessQuantile(rgSet)
mSetRaw <- preprocessRaw(rgSet)

pdf(plot2)
densityPlot(getBeta(mSetSq), main="Normalized density", legend=FALSE)
dev.off()

detP <- detP[match(featureNames(mSetSq),rownames(detP)),]
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
mSetSqFlt <- mSetSq[keep,]

cat("Filter out snps\n")
keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
mSetSqFlt <- mSetSqFlt[keep,]
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)

cat("Filter out sex chromosomes\n")
# filter out pos on sex chr
xReactiveProbes <- read.csv(file=paste(dataDirectory, "48639-non-specific-probes-Illumina450k.csv", sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
mSetSqFlt <- mSetSqFlt[keep,]
bVals <- getBeta(mSetSqFlt)

cat("Start writing to ", pathBValue, "\n")
bValsOut <- cbind(data.frame(pos = rownames(bVals)), bVals)
write_delim(x = as.data.frame(bValsOut), path = pathBValue, delim = "\t", col_names = T)
</pre>
		                     </p>
	                     </div>
	                     
	                     <div id="general_step" style="margin-top:10px;border:1px solid rgb(238,232,225);">
	                     	<B style="margin:5px 5px 5px 5px;">General_step</B>
		                     <p style="margin:15px 15px 15px 15px;">
		                     	The browse page is built based on standardized classification scheme of cancer (according to Oncotree database) and different hierarchical classification of aberrant events. To browse the mutually exclusive entries related a particular cancer types (for Papillary Thyroid Cancer), please click 'Cancers'(1 in figure below), then choose the tissue origin 'Thyroid' and corresponding cancer type 'Papillary Thyroid Cancer' (2, 3 in figure below). The related mutually exclusive entries will show on table format on the top panel, as well as mutually exclusive network on the bottom panel. Similarly, the users can browse all entries associated with an interested somatic event (Genes, Chromosomal Events OR Structural Variation Events) in a similar way.
		                     </p>
                         <p style="margin:15px 15px 15px 15px;">
		                        <B ><span id="beta2wig">beta2wig.py</span></B><br/><!--代码位置 -->
		                     	<pre class="prettyprint linenums">
import re
import sys
import glob
import multiprocessing
from collections import defaultdict


"""
this script is designed to convert 450K beta values into wig file
__author__ = likai
"""


def get_files(path_files):
    print(f"Start reading from {path_files}")
    files = glob.glob(f"{path_files}/*.beta_value.txt")
    return files


def parse_ref(path_ref):
    print(f"Start reading from {path_ref}")
    cg_pos = dict()
    with open(path_ref) as f:
        for line in f:
            line_list = line.strip().split("\t")
            cg_pos[line_list[0]] = ([line_list[1]], line_list[2])
    return cg_pos


def beta2wig(each_beta, cg_pos):
    print(f"Start parsing {each_beta}")
    beta_chr_pos = defaultdict(dict)
    beta_out = re.sub(".txt", ".wig", each_beta)
    try:
        with open(each_beta) as f:
            for line in f:
                if line.startswith("pos"):
                    continue
                line_list = line.strip().split("\t")
                if line_list[0] in cg_pos:
                    this_cg_values = cg_pos[line_list[0]]
                    beta_chr_pos[f"chr{this_cg_values[0]}"][this_cg_values[1]] = line_list[1]
                else:
                    pass
                    # print(f"no pos information found for {line_list[0]}")
        print(f"Start writing to {beta_out}")
        with open(beta_out, "w+") as out:
            for each_chr in beta_chr_pos:
                out.write(f"variableStep chrom={each_chr}\n")
                this_chr_values = beta_chr_pos[each_chr]
                for each_pos in sorted(this_chr_values, key=lambda k: int(k)):
                    out.write(each_pos + "\t" + this_chr_values[each_pos] + "\n")
    except Exception as e:
        print(e)


def multi_process_run(files, cg_pos):
    pool = multiprocessing.Pool(processes=50)
    for each_file in files:
        pool.apply_async(beta2wig, (each_file, cg_pos,))
    pool.close()
    pool.join()


def main():
    path_files = sys.argv[1]
    path_ref = sys.argv[2]
    files = get_files(path_files)
    cg_pos = parse_ref(path_ref)
    multi_process_run(files, cg_pos)


if __name__ == '__main__':
    main()
</pre>
		                     </p>
                         <p style="margin:15px 15px 15px 15px;">
		                        <B ><span id="cov2wig">cov2wig.py</span></B><br/><!--代码位置 -->
		                     	<pre class="prettyprint linenums">
import re
import sys
import glob
import gzip
import multiprocessing
from collections import defaultdict


"""
this script is designed to convert WGBS/RRBS coverage files into wig file
__author__ = likai
"""


def get_files(path_cov):
    print(f"Start reading from {path_cov}")
    cov_files = glob.glob(f"{path_cov}/*.cov.gz")
    if not cov_files:
        cov_files = glob.glob(f"{cov_files}/*.cov")
    if not cov_files:
        raise ValueError("Empty cov file dir!")
    return cov_files


def cov2wig(each_cov):
    print(f"Start parsing {each_cov}")
    total_sample_cov = defaultdict(dict)
    try:
        if each_cov.endswith(".gz"):
            f = gzip.open(each_cov)
            cov_out = re.sub(".cov.gz", ".wig.gz", each_cov)
        else:
            f = open(each_cov)
            cov_out = re.sub(".cov", ".wig.gz", each_cov)
        for line in f:
            if isinstance(line, bytes):
                line = line.decode()
            line_list = line.strip().split("\t")
            total_sample_cov[line_list[0]][line_list[1]] = line_list[3]
        print(f"Start writing to {cov_out}")
        # to save space, we output gzip format
        with gzip.open(cov_out, "wb") as out:
            for each_chr in total_sample_cov:
                out.write(f"variableStep chrom={each_chr}\n".encode())
                this_chr_values = total_sample_cov[each_chr]
                for each_pos in sorted(this_chr_values, key=lambda k: int(k)):
                    out.write(each_pos.encode() + "\t".encode() + this_chr_values[each_pos].encode() + "\n".encode())
        print(f"Finishing writing to {cov_out}")
    except Exception as e:
        print(e)


def multi_process_run(cov_files):
    pool = multiprocessing.Pool(processes=50)
    for each_file in cov_files:
        pool.apply_async(cov2wig, (each_file,))
    pool.close()
    pool.join()


def main():
    path = sys.argv[1]  # path to cov files
    cov_files = get_files(path)
    multi_process_run(cov_files)


if __name__ == '__main__':
    main()
</pre>
		                     </p>
                         <p style="margin:15px 15px 15px 15px;">
		                        <B ><span id="fastqc0">fastqc0.py</span></B><br/><!--代码位置 -->
		                     	<pre class="prettyprint linenums">
import re
import os
import glob
import argparse
import subprocess
import multiprocessing


def get_files(path_fastq):
    print(f"Start reading from {path_fastq}")
    fastq_files = glob.glob(f"{path_fastq}/*.fastq")
    return fastq_files


def run_fastqc(fastq, qc_result_dir, log_dir):
    # we run fastqc for each raw fastq sequence file
    prefix = os.path.basename(fastq).split(".")[0]
    print(f"  Fastqc: {prefix}")
    this_sample_result_dir = f"{qc_result_dir}/{prefix}"
    this_sample_log = f"{log_dir}/{prefix}.raw_qc.log"
    if not os.path.exists(this_sample_result_dir):
        os.mkdir(this_sample_result_dir)
    if fastq.endswith(".sra_1.fastq"):
        fastq_file_1 = fastq
        fastq_file_2 = f"{os.path.dirname(fastq)}/{prefix}.sra_2.fastq"
        cmd = f"fastqc -q --extract -o {this_sample_result_dir} {fastq_file_1} {fastq_file_2}"
    else:
        cmd = f"fastqc -q --extract -o {this_sample_result_dir} {fastq}"
    try:
        subprocess.check_output(cmd, shell=True)
    except Exception as e:
        print(e)
        with open(this_sample_log, "w+") as out:
            out.write(str(e) + "\n")


def multi_process_run(fastq_files, num_process, qc_result_dir, log_dir):
    pool = multiprocessing.Pool(processes=int(num_process))
    for fastq in fastq_files:
        if not fastq.endswith(".sra_2.fastq"):
            pool.apply_async(run_fastqc, (fastq, qc_result_dir, log_dir,))
    pool.close()
    pool.join()


def get_fastqc_zip_file(qc_result_dir):
    print(f"Start searching for zip results from {qc_result_dir}")
    zip_files = glob.glob(f"{qc_result_dir}/*/*.zip")
    multiqc_file_path = f"{qc_result_dir}/multiqc_sample.txt"
    with open(multiqc_file_path, "w+") as out:
        out.write("\n".join(zip_files) + "\n")
    return multiqc_file_path


def multiqc(multiqc_file_path, n, qc_result_dir):
    # we use multiqc to collect results by fastqc
    print("Start multiQC")
    cmd = f"multiqc -n {n} -s -o {qc_result_dir} -l {multiqc_file_path}"
    try:
        subprocess.check_output(cmd, shell=True)
    except Exception as e:
        print(e)


def truncate(qc_result_dir, path_out, path_fastq):
    print(f"Start searching for truncated files from {qc_result_dir}")
    samples = list(map(lambda k: re.sub(".fastq", "", os.path.basename(k)), glob.glob(f"{path_fastq}/*.fastq")))
    zip_samples = list(map(lambda k: re.sub("_fastqc.zip", "", os.path.basename(k)), glob.glob(f"{qc_result_dir}/*/*.zip")))
    truncate_files = []
    for each_sample in samples:
        if each_sample not in zip_samples:
            truncate_files.append(each_sample)
    if truncate_files:
        print(f"Start writing to {path_out}")
        with open(path_out, "w+") as out:
            out.write("\n".join(truncate_files) + "\n")


def main():
    parser = argparse.ArgumentParser(description="A program to do FastQC and MultiQC for raw fastq sequences")
    parser.add_argument("-i", action="store", dest="path_fastq")
    parser.add_argument("-p", action="store", dest="num_process", help="processes to use in the program", default=40)
    parser.add_argument("-q", action="store", dest="qc_result_dir", help="path to put qc result in")
    parser.add_argument("-l", action="store", dest="log_dir", help="path to put log files")
    parser.add_argument("-n", action="store", dest="n", help="project name to use in MultiQC output[bioproject.13]")
    results = parser.parse_args()
    if not all([results.path_fastq, results.num_process, results.qc_result_dir, results.log_dir, results.n]):
        print("too few arguments, type -h for more information!")
        exit(-1)
    path_fastq = results.path_fastq if not results.path_fastq.endswith("/") else results.path_fastq.rstrip("/")
    num_process = results.num_process
    qc_result_dir = results.qc_result_dir if not results.qc_result_dir.endswith("/") else results.qc_result_dir.rstrip("/")
    log_dir = results.log_dir if not results.log_dir.endswith("/") else results.log_dir.rstrip("/")
    n = results.n
    if not os.path.exists(qc_result_dir):
        os.makedirs(qc_result_dir)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    fastq_files = get_files(path_fastq)
    multi_process_run(fastq_files, num_process, qc_result_dir, log_dir)
    multiqc_file_path = get_fastqc_zip_file(qc_result_dir)
    multiqc(multiqc_file_path, n, qc_result_dir)
    path_truncate_out = f"{qc_result_dir}/truncate_fastq_files.txt"
    truncate(qc_result_dir, path_truncate_out, path_fastq)


if __name__ == '__main__':
    main()
</pre>
		                     </p>
                         <p style="margin:15px 15px 15px 15px;">
		                        <B ><span id="get_reads_samples_multi">get_reads_samples_multi.py</span></B><br/><!--代码位置 -->
		                     	<pre class="prettyprint linenums">
import os
import sys
import glob
import subprocess
import linecache
import multiprocessing


"""
this script is designed to count reads and samples in our project
__author__ = likai
"""


def get_samples(path):
    print(f"Start reading from {path}")
    files = glob.glob(f"{path}/*.fastq")
    if not files:
        files = glob.glob(f"{path}/*.gz")
        if not files:
            files = glob.glob(f"{path}/*/*.gz")
    return files


def count_reads(fastq_file):
    print(f"Start counting {fastq_file}")
    this_sample_tmp_out = f"{fastq_file}.tmp.txt"
    if fastq_file.endswith(".gz"):
        cmd = f"zcat {fastq_file} | wc -l > {this_sample_tmp_out}"
    else:
        cmd = f"wc -l {fastq_file} > {this_sample_tmp_out}"
    try:
        subprocess.check_output(cmd, shell=True)
    except Exception as e:
        print(e)


def run_parallel(files):
    pool = multiprocessing.Pool(processes=60)
    for fastq_file in files:
        pool.apply_async(count_reads, (fastq_file,))
    pool.close()
    pool.join()


def count_total(path, path_out, pattern):
    tmp_files = glob.glob(f"{path}/*.tmp.txt")
    if not tmp_files:
        tmp_files = glob.glob(f"{path}/*/*.tmp.txt")
    total_reads = 0
    for each_file in tmp_files:
        this_sample_reads = int(linecache.getlines(each_file)[0].strip().split()[0])
        total_reads += this_sample_reads
    for each_file in tmp_files:
        os.remove(each_file)
    samples = list(map(lambda k: os.path.basename(k).split(f"{pattern}")[0], tmp_files))
    samples = list(set(samples))
    with open(path_out, "w+") as out:
        out.write(f"total reads: {total_reads}\n")
        out.write(f"total samples: {len(samples)}")


def main():
    path = sys.argv[1]  # path to fastq files
    path_out = sys.argv[2]  # bioproject_13.reads_samples.txt
    pattern = sys.argv[3]  # pattern to split fastq names
    files = get_samples(path)
    run_parallel(files)
    count_total(path, path_out, pattern)


if __name__ == '__main__':
    main()
</pre>
		                     </p>
	                     </div>
                     	<div id="MCTA-seq" style="margin-top:10px;border:1px solid rgb(238,232,225);">
	                     	<B style="margin:5px 5px 5px 5px;">MCTA-seq</B>
												 <p style="margin:15px 15px 15px 15px;">
	                     		MutExGenome results are organized in a data table, with a mutually exclusive relationship record on each line containing “Tissue Origin”, “Cancer Type”, “Subtype”, “Gene Symbol”, “Aberrance Type”, “Method” and “RRBS_WGBS”. Users can click on the different download formats button to download the data table, and the “RRBS_WGBS” button to view detailed information of the specific entry.
	                     	</p>
		                     <p style="margin:15px 15px 15px 15px;">
		                        <B ><span id="ch3deleate">ch3deleate.pl</span></B><br/><!--代码位置 -->
		                     	<pre class="prettyprint linenums">
#!/usr/bin/perl
use strict;
use warnings;

my ($line1,$line2,$line3,$line4,$count);
my $usage = "perl $0 <IN> <OUT>
<IN>:fastq.gz
<OUT>:fastq.gz";
open IN,"zcat $ARGV[0]|" or die $!;
open OUT,"| gzip > $ARGV[1]" or die $!;
while (<IN>)
{
	$count=0;
	chomp;
	$line1 = $_;
	chomp($line2=<IN>);
	chomp($line3=<IN>);
	chomp($line4=<IN>);
	while ($line2 =~m/CA/g)
	{
		$count++;
	}
	while ($line2 =~m/CT/g)
        {
                $count++;
        }
	while ($line2 =~m/CC/g)
        {
                $count++;
        }
	next if ($count >= 3);
	print OUT "$line1\n$line2\n$line3\n$line4\n";
}
close IN;
</pre>
		                     </p>
                         <p style="margin:15px 15px 15px 15px;">
		                        <B ><span id="QCnolane">QCnolane.pl</span></B><br/><!--代码位置 -->
		                     	<pre class="prettyprint linenums">
#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
#use PerlIO::gzip;


####
#
#	Description: This sript is for Quality Control of NGS data produced via Illumia platform.
#	Main Function: 
#	(1) too much N bases. 
#	(2) trim adaptors
#	(3) delete reads with length lt 37bp after adaptor timming
#	(4) too much low quality bases
#	Author: Ping Zhu
#	Date: 2014-02-08
#	
###


#--------------------------help and options guide--------------------------#
my $usage = << USAGE;
Usage
	perl	$0
	indir	< the inout dir of the samples>
	outdir	< the output dir of the samples>
	sample	< sample>
	quality < solexa-quals/phred64-quals> [ default 33 ]
	end     < Single end = 1/Pair end = 2> [ default 2 ]
	N_rate  < N_rate> [ default 0.1 ]
	Qmin    < Qmin> [ default (quality + 5) ]
	Qrate   < Qrate> [ default 0.5 ]
USAGE

my ($indir, $outdir, $sample, $quality, $end, $N_rate, $Qmin, $Qrate, $help);
GetOptions(
    "indir=s"   => \$indir,
    "outdir=s"  => \$outdir,
    "sample=s"  => \$sample,
    #	"lane=s"=>\$lane,
    "quality=i" => \$quality,
    "end=s"     => \$end,
    "N_rate:i"  => \$N_rate,
    "Qmin:i"    => \$Qmin,
    "Qrate:i"   => \$Qrate,
    "h:s"       => \$help,
);
die $usage if $help;
#die $usage unless  $indir && $outdir && $sample && $lane;
#--------------------------help and options guide--------------------------#


#-default parameters value

$quality ||= 33;
$N_rate ||= 0.1;
$Qmin ||= ($quality + 5);
$Qrate ||= 0.5;
$end ||= 2;

#总结一下QC的默认条件：
#错误率在0.32以上的碱基达到碱基数的一半；或N碱基达到碱基数的10%的序列被去掉；同时将序列中的adaptor除掉，剩下的序列长度若低于37bp，则去掉该序列。

#-trim && create the dir if not exists

$indir = trim_slash($indir);
$outdir = trim_slash($outdir);
`mkdir -p $outdir/$sample` unless (-d "$outdir/$sample");


#-create the Log file

#open LOG,">$outdir/$sample/$lane/$sample.$lane.QC.log" or die $!;
open LOG, ">$outdir/$sample/$sample.QC.log" or die $!;
print LOG "$sample filter begin at: " . `date +%Y-%m-%d,%H:%M:%S`;
#print STDERR "$sample filter begin at: ".`date +%Y-%m-%d,%H:%M:%S`;


#-get the adapters list

my @adapter = ("GATCGGAAGAGCACA", "GATCGGAAGAGCGTC");

#-the common global variabes 

my ($total_reads, $total_bases, $remanent_reads, $remanent_bases, $without_adapter_reads, $read_length, $adapter_num, $trimAdapter_num, $remove_N_num, $low_quality_num) = (0) x 10;
my (%hash_base, %hash_quality, %hash_count, %without_adapter_bases);

my ($gc_1, $Q20_1, $Q30_1, $error_1) = (0) x 4;

#-for pair end sequencing only

my ($gc_2, $Q20_2, $Q30_2, $error_2, $remove_duplication_num) = (0) x 5;

#-get the $Q20 and $Q30

my ($Q20, $Q30);
($quality == 64) ? (($Q20, $Q30) = (84, 94)) : (($Q20, $Q30) = (53, 63));


#__________________________________________________________processing begin_______________________________________________________________#


#____________________________________for the default pair end sequencing______________________________________#

if ($end == 2) {

    #-get the input files

    chomp(my $file_1 = `ls $indir/$sample/*_R1.fastq.gz`);
    chomp(my $file_2 = `ls $indir/$sample/*_R2.fastq.gz`);

    #-open the input files

    open IN_1, "zcat $file_1 |" or die $!;
    open IN_2, "zcat $file_2 |" or die $!;

    #-open the output files

    open OUT_1, "| gzip >$outdir/$sample/$sample.1.clean.fq.gz" or die $!;
    open OUT_2, "| gzip >$outdir/$sample/$sample.2.clean.fq.gz" or die $!;


    #---------------------------read in the reads information-----------------------#

    while (1) {

        #-get the reads and corresponding information in each 4 lines

        my $line1_1 = <IN_1>;
        my $line1_2 = <IN_1>;
        my $line1_3 = <IN_1>;
        my $line1_4 = <IN_1>;

        my $line2_1 = <IN_2>;
        my $line2_2 = <IN_2>;
        my $line2_3 = <IN_2>;
        my $line2_4 = <IN_2>;

        #check the end of the file

        last unless (defined($line1_1) and defined($line2_1));

        chomp($line1_1, $line1_2, $line1_3, $line1_4, $line2_1, $line2_2, $line2_3, $line2_4);

        #count the total reads number && get the length of the first read

        $total_reads++;
        ($read_length = length($line1_2)) if ($read_length == 0);

        #-remove adapter

        my $remove_a1 = remove_adapter($line1_2, 1);
        my $remove_a2 = remove_adapter($line2_2, 2);
        if ($remove_a1 > 0 or $remove_a2 > 0) {
            if ($remove_a1 < 37 or $remove_a2 < 37) {
                $adapter_num++;
                next;
            }
            else {
                $trimAdapter_num++;
                if ($remove_a1 != 0) {
                    $line1_2 = substr($line1_2, 0, ($remove_a1 - 1));
                    $line1_4 = substr($line1_4, 0, ($remove_a1 - 1));
                }
                if ($remove_a2 != 0) {
                    $line2_2 = substr($line2_2, 0, ($remove_a2 - 1));
                    $line2_4 = substr($line2_4, 0, ($remove_a2 - 1));
                }
            }
        }

        # without_adapter_bases

        $without_adapter_bases{"1"} += length($line1_2);
        $without_adapter_bases{"2"} += length($line2_2);

        #-count bases at each site && count N content higher than %10

        my $remove_n1 = count_bases($line1_2, 0);
        my $remove_n2 = count_bases($line2_2, $read_length);
        ($remove_N_num++) if ($remove_n1 or $remove_n2);

        #-count the quality at each site && count the low quality

        my $low1 = count_quality($line1_4, \$Q20_1, \$Q30_1, 0);
        my $low2 = count_quality($line2_4, \$Q20_2, \$Q30_2, $read_length);
        ($low_quality_num++) if ($low1 or $low2);

        #-remove N content higher than %10

        if ($remove_n1 or $remove_n2) {
            next;
        }

        #-remove the low quality

        if ($low1 or $low2) {
            next;
        }

        #-count the remanent reads number

        $remanent_reads++;
        $remanent_bases += length($line1_2) + length($line2_2);

        #out put the remanent reads

        print OUT_1 "$line1_1\n$line1_2\n$line1_3\n$line1_4\n";
        print OUT_2 "$line2_1\n$line2_2\n$line2_3\n$line2_4\n";
    }

    #-caculate the reads without adapter

    $without_adapter_reads = $total_reads - $adapter_num;

    #-close the file handle

    close IN_1;
    close IN_2;
    close OUT_1;
    close OUT_2;

    #---------------------------------read in done----------------------------#


    #-----------------get the information from the variables------------------#


    #-get the total error rate && ouput the mean quality and error rate at each site

    error_rate(2);

    #-get the GC content && output the base frequency at each site

    gc_content(2);

    #-caculate a set of important rates && ouput them

    caculate_rates(2);


    #-----------------get the information from the variables done--------------#
}

#____________________________________for the default pair end sequencing done______________________________________#


#__________________________________________for the single end sequencing___________________________________________#

else {

    #-get the input files

    chomp(my $file = `ls $indir/$sample/*_R1.fastq.gz`);

    #-open the input files

    open IN_1, "<:gzip", "$file" or die $!;

    #-open the output files

    open OUT_1, ">:gzip", "$outdir/$sample/$sample.clean.fq.gz" or die $!;

    #---------------------------read in the reads information-----------------------#

    while (1) {

        #-get the reads and corresponding information in each 4 lines

        my $line1 = <IN_1>;
        my $line2 = <IN_1>;
        my $line3 = <IN_1>;
        my $line4 = <IN_1>;

        #check the end of the file

        last unless (defined($line1));
        chomp($line1, $line2, $line3, $line4);

        #count the total reads number && get the length of the first read

        $total_reads++;
        ($read_length = length($line2)) if ($read_length == 0);

        #-remove adapter

        my $remove_a1 = remove_adapter($line2, 1);
        if ($remove_a1 > 0) {
            if ($remove_a1 < 37) {
                $adapter_num++;
                next;
            }
            else {
                $trimAdapter_num++;
                $line2 = substr($line2, 0, $remove_a1 - 1);
                $line4 = substr($line4, 0, $remove_a1 - 1);
            }
        }

        # $without_adapter_bases

        $without_adapter_bases{"1"} += length($line2);

        #-count bases at each site && count N content higher than %10

        my $remove_n1 = count_bases($line2, 0);
        ($remove_N_num++) if ($remove_n1);

        #-count the quality at each site && count the low quality

        my $low1 = count_quality($line4, \$Q20_1, \$Q30_1, 0);
        ($low_quality_num++) if ($low1);

        #-remove N content higher than %10

        if ($remove_n1) {
            next;
        }

        #-remove the low quality

        if ($low1) {
            next;
        }

        #-count the remanent reads number

        $remanent_reads++;
        $remanent_bases += length($line2);

        #out put the remanent reads

        print OUT_1 "$line1\n$line2\n$line3\n$line4\n";
    }

    #-caculate the reads without adapter

    $without_adapter_reads = $total_reads - $adapter_num;

    #-close the file handle

    close IN_1;
    close OUT_1;

    #---------------------------------read in done----------------------------#


    #-----------------get the information from the variables------------------#

    #-get the total error rate && ouput the mean quality and error rate at each site

    error_rate(1);

    #-get the GC content && output the base frequency at each site

    gc_content(1);

    #-caculate a set of important rates && ouput them

    caculate_rates(1);

    #-----------------get the information from the variables done--------------#
}


#__________________________________________for the single end sequencing done_______________________________________#


#_______________________________________________plot the figures____________________________________________________#

my $X_axis;
($end == 2) ? ($X_axis = $read_length * 2) : ($X_axis = $read_length);
my $vertical_bar;
($end == 2) ? ($vertical_bar = "abline(v=$read_length,col='darkblue',lty=2)") : ($vertical_bar = "");
my $GC_figure = << FIGURE;
	gc<-read.table("$outdir/$sample/$sample.ATGC")
	site<-gc[,1]
	base_a<-gc[,4]
	base_t<-gc[,7]
	base_g<-gc[,10]
	base_c<-gc[,13]
	base_n<-gc[,16]
	total_sites<-$X_axis
	half_sites<-$read_length/2
	pdf("$outdir/$sample/$sample.ATGC.pdf",width=8,height=6)
	plot(site,base_a,xlim=c(0,total_sites),ylim=c(0,50),axes=FALSE,col="red",type="l",xlab="Position along reads",ylab="percent",main="Base percentage composition along reads",lty=1,lwd=1.5)
	lines(site,base_t,col="magenta",type="l",lty=2,lwd=1.5)
	lines(site,base_g,col="darkblue",type="l",lty=4,lwd=1.5)
	lines(site,base_c,col="green",type="l",lty=5,lwd=1.5)
	lines(site,base_n,col="cyan3",type="l",lty=6,lwd=1.5)
	legend("topright",legend=c("A","T","G","C","N"),col=c("red","magenta","darkblue","green","cyan3"),lty=c(1,2,4,5,6))
	$vertical_bar
	axis(side=1,at=seq(from=0,to=total_sites,by=half_sites))
	axis(side=2,at=seq(from=0,to=50,by=10))
	dev.off()
FIGURE
my $meanQ_errorR = << FIGURE;
	table<-read.table("$outdir/$sample/$sample.mean_quality")
        site<-table[,1]
        quality<-table[,2]
        error<-table[,3]
        total_sites<-$X_axis
        pdf("$outdir/$sample/$sample.mean_quality.pdf",width=8,height=6)
        plot(site,quality,xlim=c(0,total_sites),ylim=c(0,40),axes=FALSE,col="red",type="p",pch=".",cex=1.5,xlab="Position along reads",ylab="Quality",main="Distribution of qualities")
        axis(side=1,at=seq(from=0,to=total_sites,by=20))
        axis(side=2,at=seq(from=0,to=40,by=10))
        abline(h=20,col="darkblue",lty=2)
        abline(v=seq(0,total_sites, by=10),col="darkblue",lty=3 )

        pdf("$outdir/$sample/$sample.ErrorRate.pdf",width=8,height=6)
        plot(site,error,xlim=c(0,total_sites),col="red",type="h",xlab="Position along reads",ylab="% Error-Rate")
        axis(side=1,at=seq(from=0,to=total_sites,by=20))
        abline(v=seq(0,total_sites, by=10),col="darkblue",lty=3 )
        dev.off()
FIGURE
open R, "|R --vanilla --slave" or die $!;
print R $GC_figure;
print R $meanQ_errorR;
close R;
#_______________________________________________plot the figures done_____________________________________________#


#__________________________________________________________processing done_______________________________________________________________#


#________________________________________________Subrutines begin___________________________________________________#

#-dir trimming

sub trim_slash {
    my ($dir) = @_;
    ($dir =~ /\/$/) ? ($dir =~ s/\/$//) : ($dir = $dir);
    return $dir;
}


#-remove adapter

sub remove_adapter {
    my ($seq, $n) = @_;
    my $adapt = \@adapter;
    $n--;
    while ($seq =~ /$$adapt[$n]/g) {
        my $pos = pos($seq) - length($$adapt[$n]) + 1;
        return $pos;
    }
    return 0;
}


#-count the quality at each site && count the low quality

sub count_quality {
    my ($seq, $Q_20, $Q_30, $start_site) = @_;
    my ($i, $low_q_site, $base_quality) = (0) x 3;
    my $length = length($seq);
    while ($i < $length) {
        my $base_asc = substr($seq, $i, 1);
        $base_quality = ord($base_asc);
        $hash_quality{$i + $start_site} += $base_quality;
        $hash_count{$i + $start_site}++;
        $low_q_site++ if ($base_quality <= $Qmin);
        $$Q_20++ if ($base_quality >= $Q20);
        $$Q_30++ if ($base_quality >= $Q30);
        $i++;
    }
    ($low_q_site >= $length * $Qrate) ? (return 1) : (return 0);
}


#-count bases at each site

sub count_bases {
    my ($seq, $start_site) = @_;
    my $length = length($seq);
    my $i = 0;
    while ($i < $length) {
        my $base = substr($seq, $i, 1);
        $hash_base{$i + $start_site}{$base}++;
        $i++;
    }
    my $N_num = ($seq =~ tr/N/N/) + 0;
    (return 1) if ($N_num >= $length * $N_rate);
    return 0;
}


#-get the total error rate && ouput the mean quality and error rate at each site

sub error_rate {
    my ($end) = @_;
    open OUT_3, ">$outdir/$sample/$sample.mean_quality" or die $!;
    my @keys = sort {$a <=> $b} keys %hash_quality;
    my $minus;
    ($quality == 64) ? ($minus = 64) : ($minus = 33);
    my $i = 0;
    while ($i < @keys) {
        my $mean_quality = ($hash_quality{$keys[$i]} / $hash_count{$keys[$i]}) - $minus;
        my $index = 0 - ($mean_quality / 10);
        my $error_rate = (10 ** $index) * 100;
        if ($i < $read_length) {
            $error_1 += $error_rate;
        }
        else {
            $error_2 += $error_rate;
        }
        printf OUT_3 "%d\t%.5f\t%f\n", $keys[$i], $mean_quality, $error_rate;
        $i++;
    }
    close OUT_3;
}


#-get the GC content && output the base frequency at each site

sub gc_content {
    my ($end) = @_;
    open OUT_4, ">$outdir/$sample/$sample.ATGC" or die $!;
    my @keys = sort {$a <=> $b} keys %hash_base;
    my $i = 0;
    my @bases = qw/A T G C N/;
    while ($i < @keys) {
        print OUT_4 "$keys[$i]\t";
        my $j = 0;
        while ($j < @bases) {
            if (exists $hash_base{$keys[$i]}{$bases[$j]}) {
                my $frequency = ($hash_base{$keys[$i]}{$bases[$j]} / $hash_count{$keys[$i]}) * 100;
                printf OUT_4 "%s\t%d\t%.3f\t", $bases[$j], $hash_base{$keys[$i]}{$bases[$j]}, $frequency;
            }
            else {
                print OUT_4 "$bases[$j]\t0\t0\t";
            }
            $j++;
        }
        print OUT_4 "\n";

        my ($g, $c) = (0) x 2;
        ($g = $hash_base{$keys[$i]}{"G"}) if (exists($hash_base{$keys[$i]}{"G"}));
        ($c = $hash_base{$keys[$i]}{"C"}) if (exists($hash_base{$keys[$i]}{"C"}));
        ($i < $read_length) ? ($gc_1 += $g + $c) : ($gc_2 += $g + $c);
        $i++;
    }
    close OUT_4;
}


#-caculate a set of important rates && ouput them

sub caculate_rates {
    my ($end) = @_;
    my ($gc_rate_2, $Q20_rate_2, $Q30_rate_2, $error_rate_2, $duplication_rate);
    $total_bases = $total_reads * $read_length;
    my $gc_rate_1 = ($gc_1 / $without_adapter_bases{"1"}) * 100;
    my $Q20_rate_1 = ($Q20_1 / $without_adapter_bases{"1"}) * 100;
    my $Q30_rate_1 = ($Q30_1 / $without_adapter_bases{"1"}) * 100;
    my $error_rate_1 = $error_1 / 100;
    if ($end == 2) {
        $total_reads = $total_reads * 2;
        $total_bases = $total_reads * $read_length;
        $remanent_reads = $remanent_reads * 2;
        $gc_rate_2 = ($gc_2 / $without_adapter_bases{"2"}) * 100;
        $Q20_rate_2 = ($Q20_2 / $without_adapter_bases{"2"}) * 100;
        $Q30_rate_2 = ($Q30_2 / $without_adapter_bases{"2"}) * 100;
        $error_rate_2 = $error_2 / 100;
    }
    my $remanent_base_rate = $remanent_bases / $total_bases * 100;

    my $title = "Raw reads\tRaw bases\tClean reads\tClean bases\tCleanRate\tErrorRate\tQ20\tQ30\tGC content\n";
    printf LOG $title;

    my $output1 = "$total_reads\t$total_bases\t$remanent_reads\t$remanent_bases\t$remanent_base_rate\t";
    printf LOG $output1;

    if ($end == 2) {
        printf LOG "%.2f;%.2f\t%.2f;%.2f\t%.2f;%.2f\t%.2f;%.2f\n", $error_rate_1, $error_rate_2, $Q20_rate_1, $Q20_rate_2, $Q30_rate_1, $Q30_rate_2, $gc_rate_1, $gc_rate_2;
    }
    else {
        printf LOG "%.2f\t%.2f\t%.2f\t%.2f\n", $error_rate_1, $Q20_rate_1, $Q30_rate_1, $gc_rate_1;
    }
    my $output2 = "N remove $remove_N_num\nQuality remove $low_quality_num\nAdapter remove $adapter_num\tTrimed Adapter $trimAdapter_num\n";
    printf LOG $output2;

    printf LOG "$sample filter end at: " . `date +%Y-%m-%d,%H:%M:%S`;
    close LOG;
}

#________________________________________________Subrutines done___________________________________________________#
</pre>
		                     </p>
                         <p style="margin:15px 15px 15px 15px;">
		                        <B ><span id="rm_firstx_leny">rm_firstx_leny.pl</span></B><br/><!--代码位置 -->
		                     	<pre class="prettyprint linenums">
#! /usrt/bin/perl -w
use strict;
use warnings FATAL => 'all';


# remove forward xbp, backward zbp, remain more than ybp
my $usage = "perl $0 <IN1> <OUT> x y z";
die $usage unless @ARGV==5;
open IN,"zcat $ARGV[0] |" or die $!;
open OUT,"| gzip > $ARGV[1]" or die $!;
while(<IN>)
{
	chomp;
	my $line1=$_;
	chomp(my $line2=<IN>);
	chomp(my $line3=<IN>);
	chomp(my $line4=<IN>);
	my $len_2=length($line2);
	if($len_2 >=$ARGV[2]+$ARGV[4]+$ARGV[3])
	{
		my $second=substr($line2,$ARGV[2],$len_2-$ARGV[2]-$ARGV[4]);
		my $line4_1=substr($line4,$ARGV[2],$len_2-$ARGV[2]-$ARGV[4]);
		my $len=length($second);
		#my $third=substr($second,0,$len_2-24);
		#my $line4_2=substr($line4_1,0,$len_2-24);
		print OUT "$line1\n$second\n$line3\n$line4_1\n";
	}
	else
	{
		next;
	}
}
close IN;
close OUT;
</pre>
		                     </p>
                         <p style="margin:15px 15px 15px 15px;">
		                        <B ><span id="trim_mapping_bio_19">trim_mapping_bio_19.py</span></B><br/><!--代码位置 -->
		                     	<pre class="prettyprint linenums">
import os
import glob
import argparse
import subprocess
import multiprocessing


"""
this script is written to align fastq files by MCTA-seq. Before aligning two reads separately, 
we use three custom perl scripts to trim fastq files.
__author__ = likai 
"""


def get_files(path_fastq):
    """
    get all the fastq files from MCTA-seq, these files need to be ended with .gz
    :param path_fastq:
    :return: python list with all gzipped fastq files
    """
    print(f"Start reading from {path_fastq}")
    fastq_files = glob.glob(f"{path_fastq}/*/*.fastq.gz")
    return fastq_files


def trim_map_file(fastq, path_fastq, clean_dir, log_dir, mapping_dir, bismark_index, bowtie1_path, mapping_result_dir):
    """
    we first use QCnolane.pl to trim both fastq files, then use rm_firstx_leny.pl
    to remove forward xbp, backward zbp, remain more than ybp in the trimmed fastq
    files. Bismark is used to align read1 and read2 separately.
    :param fastq: path to the input gzipped fastq file
    :param path_fastq: directory to all fastq files
    :param clean_dir: directory to put trimmed fastq files
    :param log_dir: directory to put log files
    :param mapping_dir: directory to put aligned bam files
    :param bismark_index: path to bismark index files
    :param bowtie1_path: path to bowtie1
    :param mapping_result_dir: directory to put bismark log files
    :return: none
    """
    prefix = os.path.basename(fastq).split("_")[0]
    this_sample_error_log = f"{log_dir}/{prefix}.error.log"
    cmd_trim = f"perl QCnolane.pl --indir {path_fastq} --outdir {clean_dir} --sample {prefix}"
    clean_fq_1 = f"{clean_dir}/{prefix}/{prefix}.1.clean.fq.gz"
    clean_fq_2 = f"{clean_dir}/{prefix}/{prefix}.2.clean.fq.gz"
    filter_fq_1 = f"{clean_dir}/{prefix}/{prefix}.1.filter.fq.gz"
    filter_fq_2 = f"{clean_dir}/{prefix}/{prefix}.2.filter.fq.gz"
    nonCG_fq_2 = f"{clean_dir}/{prefix}/{prefix}.2.nonCG.fq.gz"
    bismark_output_1 = f"{mapping_dir}/{prefix}_1"
    bismark_output_2 = f"{mapping_dir}/{prefix}_2"
    if not os.path.exists(bismark_output_1):
        os.makedirs(bismark_output_1)
    if not os.path.exists(bismark_output_2):
        os.makedirs(bismark_output_2)
    cmd_filter_2 = f"perl rm_firstx_leny.pl {clean_fq_2} {filter_fq_2} 6 20 12"
    cmd_filter_1 = f"perl rm_firstx_leny.pl {clean_fq_1} {filter_fq_1} 12 20 6"
    bismark_log_1 = f"{mapping_result_dir}/{prefix}/{prefix}.1.mapping.txt"
    bismark_log_2 = f"{mapping_result_dir}/{prefix}/{prefix}.2.mapping.txt"
    if not os.path.exists(f"{mapping_result_dir}/{prefix}"):
        os.makedirs(f"{mapping_result_dir}/{prefix}")
    cmd_ch3 = f"perl ch3deleate.pl {filter_fq_2} {nonCG_fq_2}"
    # we use bismark to align each read separately
    cmd_bismark_1 = f"bismark --bowtie1 --non_directional --fastq --phred33-quals --temp_dir {log_dir} " \
                    f"--path_to_bowtie {bowtie1_path} --output_dir {bismark_output_1} {bismark_index}" \
                    f" {filter_fq_1} > {bismark_log_1} 2>&1"
    cmd_bismark_2 = f"bismark --bowtie1 --non_directional --fastq --phred33-quals --temp_dir {log_dir} " \
                    f"--path_to_bowtie {bowtie1_path} --output_dir {bismark_output_2} {bismark_index} " \
                    f"{nonCG_fq_2} > {bismark_log_2} 2>&1"
    bam1 = f"{bismark_output_1}/{prefix}.1.filter_bismark.bam"
    bam2 = f"{bismark_output_2}/{prefix}.2.nonCG_bismark.bam"
    cmd_samtools_sort_1 = f"samtools sort -@ 20 {bam1} -o {os.path.dirname(bam1)}/{prefix}.1.sort.bam"
    cmd_samtools_sort_2 = f"samtools sort -@ 20 {bam2} -o {os.path.dirname(bam2)}/{prefix}.2.sort.bam"
    cmd_index_1 = f"samtools index {os.path.dirname(bam1)}/{prefix}.1.sort.bam"
    cmd_index_2 = f"samtools index {os.path.dirname(bam2)}/{prefix}.2.sort.bam"
    try:
        print(f"    Trimming:  {prefix}")
        subprocess.check_output(cmd_trim, shell=True)
        print(f"    Filter fq 2:  {prefix}")
        subprocess.check_output(cmd_filter_2, shell=True)
        print(f"    Filter fq 1:  {prefix}")
        subprocess.check_output(cmd_filter_1, shell=True)
        print(f"    Filter out CG of fq 2:  {prefix}")
        subprocess.check_output(cmd_ch3, shell=True)
        print(f"    Aligning fq 1:  {prefix}")
        subprocess.check_output(cmd_bismark_1, shell=True)
        print(f"    Aligning fq 2:  {prefix}")
        subprocess.check_output(cmd_bismark_2, shell=True)
        print(f"    Sorting bam1:  {prefix}")
        subprocess.check_output(cmd_samtools_sort_1, shell=True)
        print(f"    Sorting bam2:  {prefix}")
        subprocess.check_output(cmd_samtools_sort_2, shell=True)
        print(f"    Indexing bam1:  {prefix}")
        subprocess.check_output(cmd_index_1, shell=True)
        print(f"    Indexing bam2:  {prefix}")
        subprocess.check_output(cmd_index_2, shell=True)
    except Exception as e:
        print(e)
        with open(this_sample_error_log, "a") as out:
            out.write(str(e) + "\n")


def multi_process_run(fastq_files, num_process, path_fastq, clean_dir, log_dir, mapping_dir, bismark_index,
                      bowtie1_path, mapping_result_dir):
    pool = multiprocessing.Pool(processes=int(num_process))
    print("Start checking files")
    for each_fastq in fastq_files:
        prefix = os.path.basename(each_fastq).split("_")[0]
        if not os.path.exists(f"{os.path.dirname(each_fastq)}/{prefix}_R1.fastq.gz") or not os.path.exists(
                f"{os.path.dirname(each_fastq)}/{prefix}_R2.fastq.gz"):
            print(f"{prefix} not completed!")
            exit(-1)
    for fastq in fastq_files:
        if not fastq.endswith("_R2.fastq.gz"):
            pool.apply_async(trim_map_file, (
                fastq, path_fastq, clean_dir, log_dir, mapping_dir, bismark_index, bowtie1_path, mapping_result_dir))
    pool.close()
    pool.join()


def main():
    parser = argparse.ArgumentParser(description="A program to map raw fastq for MCTA-seq")
    parser.add_argument("-i", action="store", dest="path_fastq")
    parser.add_argument("-p", action="store", dest="num_process", help="processes to use in the program", default=40)
    parser.add_argument("-m", action="store", dest="mapping_result_dir", help="path to put mapping result in")
    parser.add_argument("-l", action="store", dest="log_dir", help="path to put log files")
    parser.add_argument("-b", action="store", dest="mapping_dir", help="path to put bam files")
    parser.add_argument("-t", action="store", dest="bowtie1_path", help="path to bowtie1 executable")
    parser.add_argument("-d", action="store", dest="path_index", help="path to bismark/bowtie1 indexes")
    parser.add_argument("-c", action="store", dest="clean_dir", help="path to put trimmed fastq files")
    results = parser.parse_args()
    if not all([results.path_fastq, results.num_process, results.mapping_result_dir, results.log_dir,
                results.mapping_dir, results.bowtie1_path, results.path_index, results.clean_dir]):
        print("too few arguments, type -h for more information!")
        exit(-1)
    path_fastq = results.path_fastq if not results.path_fastq.endswith("/") else results.path_fastq.rstrip("/")
    num_process = results.num_process
    mapping_result_dir = results.mapping_result_dir if not results.mapping_result_dir.endswith(
        "/") else results.mapping_result_dir.rstrip("/")
    log_dir = results.log_dir if not results.log_dir.endswith("/") else results.log_dir.rstrip("/")
    mapping_dir = results.mapping_dir if not results.mapping_dir.endswith("/") else results.mapping_dir.rstrip("/")
    bowtie1_path = results.bowtie1_path if not results.bowtie1_path.endswith("/") else results.bowtie1_path.rstrip("/")
    bismark_index = results.path_index if not results.path_index.endswith("/") else results.path_index.rstrip("/")
    clean_dir = results.clean_dir if not results.clean_dir.endswith("/") else results.clean_dir.rstrip("/")
    for each_dir in [mapping_result_dir, log_dir, mapping_dir, clean_dir]:
        if not os.path.exists(each_dir):
            os.makedirs(each_dir)
    fastq_files = get_files(path_fastq)
    multi_process_run(fastq_files, num_process, path_fastq, clean_dir, log_dir,
                      mapping_dir, bismark_index, bowtie1_path, mapping_result_dir)


if __name__ == '__main__':
    main()
</pre>
		                     </p>
	                     </div>
	                     <div id="Nucelosome" style="margin-top:10px;border:1px solid rgb(238,232,225);">
	                     	<B style="margin:5px 5px 5px 5px;">Nucelosome</B>
	                     	<p style="margin:15px 15px 15px 15px;">
	                     		MutExGenome results are organized in a data table, with a mutually exclusive relationship record on each line containing “Tissue Origin”, “Cancer Type”, “Subtype”, “Gene Symbol”, “Aberrance Type”, “Method” and “RRBS_WGBS”. Users can click on the different download formats button to download the data table, and the “RRBS_WGBS” button to view detailed information of the specific entry.
	                     	</p>
                         <p style="margin:15px 15px 15px 15px;">
		                        <B ><span id="atacCallPeak">atacCallPeak.sh</span></B><br/><!--代码位置 -->
		                     	<pre class="prettyprint linenums">

</pre>
		                     </p>
                         <p style="margin:15px 15px 15px 15px;">
		                        <B ><span id="atacMapTrim">atacMapTrim.sh</span></B><br/><!--代码位置 -->
		                     	<pre class="prettyprint linenums">

</pre>
		                     </p>
	                     </div>
	                     <div id="RRBS_WGBS" style="margin-top:10px;border:1px solid rgb(238,232,225);">
	                     	<B style="margin:5px 5px 5px 5px;">RRBS_WGBS</B>
	                     	<p style="margin:15px 15px 15px 15px;">Users can obtain detailed information on each mutually exclusive entry by clicking on the “RRBS_WGBS” button of the data table, which opens a page containing comprehensive information of the single record and provides links to external annotation web sites, as well as mutually exclusive heat map if possible.</p>
	                     	<p style="margin:15px 15px 15px 15px;">
		                        <B ><span id="extract_cov_li_20">extract_cov_li_20.py</span></B><br/><!--代码位置 -->
		                     	<pre class="prettyprint linenums">
import os
import glob
import argparse
import subprocess
import multiprocessing


def get_bam_files(mapping_dir):
    print(f"Start reading from {mapping_dir}")
    bam_files = glob.glob(f"{mapping_dir}/*/*.bam")
    return bam_files


def bam2cov(bam, cov_dir, log_dir):
    prefix = os.path.basename(bam).split(".")[0]
    this_cov_dir = f"{cov_dir}/{prefix}"
    this_log_dir = f"{log_dir}/{prefix}.bismark_methylation_extractor.log"
    this_error_log = f"{log_dir}/{prefix}.bismark_methylation_extractor.error.log"
    if not os.path.exists(this_cov_dir):
        os.makedirs(this_cov_dir)
    if bam.endswith("_bt2.bam"):
        cmd_cov = f"bismark_methylation_extractor -s {bam} --bedGraph --counts -o {this_cov_dir} > {this_log_dir} 2>&1"
    elif bam.endswith("_bt2_pe.bam"):
        cmd_cov = f"bismark_methylation_extractor -p {bam} --bedGraph --counts -o {this_cov_dir} > {this_log_dir} 2>&1"
    else:
        cmd_cov = "echo -e '\tInvalid bam, neither paired nor single!'"
    try:
        print(f"    Extracting coverage: {prefix}")
        subprocess.check_output(cmd_cov, shell=True)
    except Exception as e:
        print(e)
        with open(this_error_log, "w+") as out:
            out.write(str(e) + "\n")


def multi_process_run(bam_files, num_process, cov_dir, log_dir):
    pool = multiprocessing.Pool(processes=int(num_process))
    for bam in bam_files:
        pool.apply_async(bam2cov, (bam, cov_dir, log_dir,))
    pool.close()
    pool.join()


def main():
    parser = argparse.ArgumentParser(description="A program to map raw fastq and remove duplicates")
    parser.add_argument("-i", action="store", dest="mapping_dir")
    parser.add_argument("-p", action="store", dest="num_process", help="processes to use in the program", default=10)
    parser.add_argument("-c", action="store", dest="cov_dir", help="path to put mapping result in")
    parser.add_argument("-l", action="store", dest="log_dir", help="path to put log files")
    results = parser.parse_args()
    if not all([results.mapping_dir, results.num_process, results.cov_dir, results.log_dir]):
        print("too few arguments, type -h for more information")
        exit(-1)
    mapping_dir = results.mapping_dir if not results.mapping_dir.endswith("/") else results.mapping_dir.rstrip("/")
    num_process = results.num_process
    cov_dir = results.cov_dir if not results.cov_dir.endswith("/") else results.cov_dir.rstrip("/")
    log_dir = results.log_dir if not results.log_dir.endswith("/") else results.log_dir.rstrip("/")
    for each_dir in [cov_dir, log_dir]:
        if not os.path.exists(each_dir):
            os.makedirs(each_dir)
    bam_files = get_bam_files(mapping_dir)
    multi_process_run(bam_files, num_process, cov_dir, log_dir)


if __name__ == '__main__':
    main()
</pre>
		                     </p>
                         <p style="margin:15px 15px 15px 15px;">
		                        <B ><span id="get_meta_sql_li_20">get_meta_sql_li_20.py</span></B><br/><!--代码位置 -->
		                     	<pre class="prettyprint linenums">
import sys
from collections import defaultdict


def get_qc_result(path_qc_result):
    print(f"Start reading frm {path_qc_result}")
    sample_qc_results = defaultdict(dict)
    with open(path_qc_result) as f:
        for line in f:
            if line.startswith("Sample"):
                continue
            line_list = line.strip().split("\t")
            sample_srr = line_list[0].split(".")[0]
            sample_qc_results[sample_srr]["adapter_content"] = line_list[1]
            sample_qc_results[sample_srr]["duplication_level"] = line_list[3]
            sample_qc_results[sample_srr]["base_quality"] = line_list[6]
            sample_qc_results[sample_srr]["basic_statistics"] = line_list[10]
            sample_qc_results[sample_srr]["overrepresented_sequences"] = line_list[15]
    return sample_qc_results


def get_mapping_results(path_mapping_result, sample_qc_results):
    print(f"Start reading from {path_mapping_result}")
    with open(path_mapping_result) as f:
        for line in f:
            if line.startswith("Sample"):
                continue
            line_list = line.strip().split("\t")
            sample_srr = line_list[0].split(".")[0]
            if sample_srr in sample_qc_results:
                sample_qc_results[sample_srr]["alignment_rate"] = line_list[1]
            else:
                raise ValueError(f"sample {sample_srr} has no mapping results!")
    return sample_qc_results


def parse_meta(path_meta, sample_qc_results):
    print(f"Start parsing {path_meta}")
    with open(path_meta) as f:
        for line in f:
            if line.startswith("sample"):
                continue
            line_list = line.strip().split("\t")
            sample_srr = line_list[0]
            if sample_srr in sample_qc_results:
                sample_qc_results[sample_srr]["GSM"] = line_list[0]
                sample_qc_results[sample_srr]["disease_state"] = line_list[1]
                sample_qc_results[sample_srr]["PMID"] = line_list[2]
                sample_qc_results[sample_srr]["method"] = line_list[3]
                sample_qc_results[sample_srr]["sample_source"] = line_list[5]
                sample_qc_results[sample_srr]["term"] = line_list[9]
                sample_qc_results[sample_srr]["gender"] = line_list[10]
                sample_qc_results[sample_srr]["age"] = line_list[11]
                sample_qc_results[sample_srr]["stage"] = line_list[12]
            else:
                print(f"sample {sample_srr} is not in liquid_20!!")
    return sample_qc_results


def write(sample_qc_results, path_out):
    print(f"Start writing to {path_out}")
    titles = ["GSM", "disease_state", "PMID", "method", "sample_source", "term",
              "gender", "age", "stage", "adapter_content", "duplication_level",
              "base_quality", "basic_statistics", "overrepresented_sequences", "alignment_rate"]
    with open(path_out, "w+") as out:
        out.write("sample\t" + "\t".join(titles) + "\n")
        for each_sample in sample_qc_results:
            this_sample_values = [sample_qc_results[each_sample][each_item] for each_item in titles]
            out.write(each_sample + "\t" + "\t".join(this_sample_values) + "\n")


def main():
    path_qc_result = sys.argv[1]  # multiqc_fastqc.txt
    path_mapping_result = sys.argv[2]  # multiqc_bowtie2.txt
    path_meta = sys.argv[3]  # cfDNA_samples_meta_plasma.xls
    path_out = sys.argv[4]  # bio_project_13_meta_qc_to_sql.xls
    sample_qc_results = get_qc_result(path_qc_result)
    sample_qc_results = get_mapping_results(path_mapping_result, sample_qc_results)
    sample_qc_results = parse_meta(path_meta, sample_qc_results)
    write(sample_qc_results, path_out)


if __name__ == '__main__':
    main()
</pre>
		                     </p>
                         <p style="margin:15px 15px 15px 15px;">
		                        <B ><span id="getCG-1.0.0">getCG-1.0.0.pl</span></B><br/><!--代码位置 -->
		                     	<pre class="prettyprint linenums">
#!/usr/bin/perl
use warnings;
use strict;

# this script is written to extract CpG positions from whole genome fasta file
# the input of this script is a whole genome fastq file, foe example, hg19.fa
# you can simply run the script like perl getCG-1.0.0.pl hg19.fa CpG.hg19.txt
# __author__ == wangxinyu

my $marker = "CG";										#marker

die "Thin program help you get whole genome CpGs position.
(perl) getCG genome.fa geneome_CpG.txt
(perl) getCG hg19/ genome_CpG.txt
	(require .fa files in hg19/)\n" unless(@ARGV);
my ( $infile , $outfile ) = @ARGV;
my ($in_fh,$out_fh);

$marker = uc($marker);									#To capital letter.

if(-e $infile and -d $infile){
	print "Input is a directory.\n";
	
	open $out_fh,">",$outfile or die "$outfile open ERROR!\n";
	
	my @files = glob("$infile/*.fa");
	foreach my $chr_file(@files){
		open $in_fh,"<",$chr_file or die "$chr_file open ERROR!\n";
		my ($chr , $chr_string) = ( "" , "" );
		
		#new chrom start	
		$chr = <$in_fh>;
		chomp $chr;
		$chr =~ s/>//;
		print $chr." is running...\n";
		
		while(my $line = <$in_fh>){
			chomp $line;
			$chr_string .= uc($line);
		}
		print_marker_position($chr,$chr_string,$marker,$out_fh);
	}
}
elsif(-e $infile){
	print "Input is a file.\n";
	
	open $in_fh,"<",$infile or die "$infile open ERROR!\n";
	open $out_fh,">",$outfile or die "$outfile open ERROR!\n";
	
	my ($chr , $chr_string) = ( "" , "" );
	while(my $line = <$in_fh>){
		chomp $line;
		if($line =~ />/){
			$line =~ s/>//;
			print $line."  is running...\n";
			if($chr){
				print_marker_position($chr,$chr_string,$marker,$out_fh);
			}
			
			#new chrom start
			$chr_string = "";
			$chr = $line;
		}
		else{
			$chr_string .= uc($line);
		}
	}
	#the last chrom
	print_marker_position($chr,$chr_string,$marker,$out_fh);
	
	close $in_fh;
	close $out_fh;
}
else{
	print "INPUT ERROR!\nCheck you input:\n$infile\n";
}

#print_marker_position($chr,$chr_string,$marker,$out_fh)
sub print_marker_position{
	my ($chr,$chr_string,$marker,$out_fh) = @_;
	
	my $position = index($chr_string,$marker) + 1;
	while($position ne 0){
		print $out_fh "$chr\t$position\n";
		$position = index($chr_string,$marker,$position) + 1;
	}
}
</pre>
		                     </p>
                         <p style="margin:15px 15px 15px 15px;">
		                        <B ><span id="grep_cpg_li_20">grep_cpg_li_20.py</span></B><br/><!--代码位置 -->
		                     	<pre class="prettyprint linenums">
import glob
import gzip
import argparse
import multiprocessing


"""
this script is designed to merge CpG positions in the positive and
negative strand
__author__ = likai 
"""


def get_cov_files(cov_dir):
    """
    get all coverage file generated by bismark
    :param cov_dir: directory to coverage file by bismark
    :return: none
    """
    print(f"Start reading from {cov_dir}")
    cov_files = glob.glob(f"{cov_dir}/*/*.cov.gz")
    return cov_files


def parse_ref(path_ref):
    """
     read in the reference CpG file into a python dict
    :param path_ref: path to reference CpG position in the whole genome fasta file
    :return: python dict, key stands for chromosome + position
    """
    print("Start reading from {}".format(path_ref))
    ref_dict = dict()
    with gzip.open(path_ref) as f:
        for line in f:
            line = line.decode()
            line_list = line.strip().split()
            key = line_list[0] + "-" + line_list[-1]
            ref_dict[key] = 1
    return ref_dict


def parse_cov(path_cov, ref_dict):
    """
    merge CpG position in the positive and negative strand
    :param path_cov:
    :param ref_dict:
    :return:
    """
    try:
        print("Start reading from {}".format(path_cov))
        path_out = path_cov[:-6] + "met.cov"
        path_temp = path_cov[:-6] + "tmp.cov"
        out = open(path_out, "w+")
        temp = open(path_temp, "w+")
        out_lines = dict()
        with gzip.open(path_cov) as f:
            for line in f:
                line = line.decode()
                line_list = line.strip().split()
                key_positive = line_list[0] + "-" + line_list[1]
                key_negative = line_list[0] + "-" + str(int(line_list[1]) - 1)
                if key_positive in ref_dict:
                    out_lines[key_positive] = line_list[:3] + line_list[4:]
                elif key_negative in ref_dict:
                    if key_negative in out_lines:
                        negative_value = out_lines[key_negative]
                        negative_read = list(map(int, negative_value[3:]))
                        line_list_read = list(map(int, line_list[4:]))
                        out_read = [negative_read[0] + line_list_read[0], negative_read[1] + line_list_read[1]]
                        out_read = list(map(str, out_read))
                        out_lines[key_negative] = negative_value[:3] + out_read
                    else:
                        out_lines[key_negative] = [line_list[0], str(int(line_list[1]) - 1),
                                                   str(int(line_list[2]) - 1)] + line_list[4:]
                else:
                    temp.write(line)
        print("Start writing to {}".format(path_out))
        for each_key in out_lines:
            value = out_lines[each_key]
            out.write("\t".join(value) + "\n")
    except Exception as e:
        print(e)


def multiprocess_run(cov_files, ref_dict, num_process):
    pool = multiprocessing.Pool(processes=int(num_process))
    for each_cov in cov_files:
        pool.apply_async(parse_cov, (each_cov, ref_dict,))
    pool.close()
    pool.join()


def main():
    parser = argparse.ArgumentParser(description="A program to extract CPG from bismark coverage files")
    parser.add_argument("-i", action="store", dest="cov_dir")
    parser.add_argument("-p", action="store", dest="num_process", help="processes to use in the program", default=10)
    parser.add_argument("-f", action="store", dest="path_ref", help="path to hg19 whole genome CPGs")
    results = parser.parse_args()
    if not all([results.cov_dir, results.num_process, results.path_ref]):
        print("too few arguments, type -h for more information")
        exit(-1)
    cov_dir = results.cov_dir if not results.cov_dir.endswith("/") else results.cov_dir.rstrip("/")
    num_process = results.num_process
    path_ref = results.path_ref
    cov_files = get_cov_files(cov_dir)
    ref_dict = parse_ref(path_ref)
    multiprocess_run(cov_files, ref_dict, num_process)


if __name__ == '__main__':
    main()
</pre>
		                     </p>
                         <p style="margin:15px 15px 15px 15px;">
		                        <B ><span id="mergeFq">mergeFq.py</span></B><br/><!--代码位置 -->
		                     	<pre class="prettyprint linenums">
import os
import re
import sys
import subprocess
from itertools import chain
from collections import defaultdict


"""
this script is written to merge different runs from same samples
"""


def get_run_files(path_meta, pmid):
    """
    
    :param path_meta:
    :param pmid:
    :return:
    """
    print(f"Start reading from {path_meta}")
    sample_run = dict()
    with open(path_meta) as f:
        for line in f:
            line_list = line.strip().split("\t")
            if line_list[2] == pmid:
                urls = line_list[7]
                urls_list = list(map(lambda k: re.sub(".sra", "", os.path.basename(k)), urls.strip().split(",")))
                sample_run[line_list[0]] = urls_list
    return sample_run


def check_run_file(sample_run, path_fastq):
    """
    collect run information for each sample
    :param sample_run: run information for each sample
    :param path_fastq: directory to each fastq file
    :return:
    """
    sample_run_out = defaultdict(lambda: defaultdict(list))
    print("Start checking run files")
    for each_sample in sample_run:
        print(f" Start parsing {each_sample}")
        this_sample_runs = list(map(lambda k: [f"{path_fastq}/{k}.sra.fastq",
                                               f"{path_fastq}/{k}.sra_1.fastq",
                                               f"{path_fastq}/{k}.sra_2.fastq"], sample_run[each_sample]))
        this_sample_runs = list(chain(*this_sample_runs))
        this_sample_fastq = [each_fastq for each_fastq in this_sample_runs if os.path.exists(each_fastq)]
        if not this_sample_fastq:
            raise ValueError(f"{each_sample} has no fastq files!")
        for each_fastq in this_sample_fastq:
            if each_fastq.endswith(".sra.fastq"):
                sample_run_out[each_sample]["read"].append(each_fastq)
            elif each_fastq.endswith(".sra_1.fastq"):
                sample_run_out[each_sample]["read1"].append(each_fastq)
            elif each_fastq.endswith(".sra_2.fastq"):
                sample_run_out[each_sample]["read2"].append(each_fastq)
            else:
                raise ValueError(f"sample {each_sample} has invalid fastq file types!")
    return sample_run_out


def get_cmd_paired(sample_run_out, path_out, each_sample):
    """
    shell command to merge the paired fastq files
    :param sample_run_out:
    :param path_out:
    :param each_sample:
    :return:
    """
    this_sample_read1_files = sample_run_out[each_sample]["read1"]
    this_sample_read2_files = sample_run_out[each_sample]["read2"]
    this_sample_read1_out = f"{path_out}/{each_sample}.sra_1.fastq"
    this_sample_read2_out = f"{path_out}/{each_sample}.sra_2.fastq"
    cmd_read1 = f"cat {' '.join(this_sample_read1_files)} > {this_sample_read1_out}"
    cmd_read2 = f"cat {' '.join(this_sample_read2_files)} > {this_sample_read2_out}"
    return cmd_read1, cmd_read2


def get_cmd_single(sample_run_out, path_out, each_sample):
    """
    shell command to merge the single end fastq files
    :param sample_run_out:
    :param path_out:
    :param each_sample:
    :return:
    """
    this_sample_read_files = sample_run_out[each_sample]["read"]
    this_sample_read_out = f"{path_out}/{each_sample}.sra.fastq"
    cmd_read = f"cat {' '.join(this_sample_read_files)} > {this_sample_read_out}"
    return cmd_read


def merge_fastq(sample_run_out, path_out, log_file):
    """
    rum the merge shell command to merge different run files in each sample
    :param sample_run_out:
    :param path_out:
    :param log_file:
    :return:
    """
    try:
        log = open(log_file, "w+")
        for each_sample in sample_run_out:
            print(f"  Start merging sample {each_sample}")
            this_sample_read_types = list(sample_run_out[each_sample].keys())
            if len(this_sample_read_types) == 3:
                print(f"   sample {each_sample} has both single end and paired end files!")
                log.write(f"{each_sample}")
                read_fastq_sizes = list(map(lambda k: int(os.path.getsize(k)), sample_run_out[each_sample]["read"]))
                read1_fastq_sizes = list(map(lambda k: int(os.path.getsize(k)), sample_run_out[each_sample]["read1"]))
                read2_fastq_sizes = list(map(lambda k: int(os.path.getsize(k)), sample_run_out[each_sample]["read2"]))
                if sum(read1_fastq_sizes) + sum(read2_fastq_sizes) >= sum(read_fastq_sizes):
                    cmd_read1, cmd_read2 = get_cmd_paired(sample_run_out, path_out, each_sample)
                    print("  Using paired fastq files instead!")
                    subprocess.check_output(cmd_read1, shell=True)
                    subprocess.check_output(cmd_read2, shell=True)
                else:
                    cmd_read = get_cmd_single(sample_run_out, path_out, each_sample)
                    print("  Using single fastq files instead!")
                    subprocess.check_output(cmd_read, shell=True)
            if this_sample_read_types == ["read1", "read2"]:
                cmd_read1, cmd_read2 = get_cmd_paired(sample_run_out, path_out, each_sample)
                print(cmd_read1)
                subprocess.check_output(cmd_read1, shell=True)
                print(cmd_read2)
                subprocess.check_output(cmd_read2, shell=True)
            elif this_sample_read_types == ["read"]:
                cmd_read = get_cmd_single(sample_run_out, path_out, each_sample)
                print(cmd_read)
                subprocess.check_output(cmd_read, shell=True)
            else:
                print(f"Invalid fastq file types {each_sample}")
                log.write(f"{each_sample} invalid fastq types")
        log.close()
    except Exception as e:
        print(e)


def main():
    path_meta = sys.argv[1]  # cfDNA_samples_meta_plasma.xls
    path_fastq = sys.argv[2]  # /share/pub/lik/cfDNA_database/projects/liquid.20/fastq
    path_out = sys.argv[3]  # /share/pub/lik/cfDNA_database/projects/liquid.20/fastq_merge
    pmid = sys.argv[4]  # 28263317 PMID number
    log_file = sys.argv[5]  # /share/pub/lik/cfDNA_database/projects/liquid.20/log/invalid_fastq_sample.txt
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    sample_run = get_run_files(path_meta, pmid)
    sample_run_out = check_run_file(sample_run, path_fastq)
    merge_fastq(sample_run_out, path_out, log_file)


if __name__ == '__main__':
    main()
</pre>
		                     </p>
                         <p style="margin:15px 15px 15px 15px;">
		                        <B ><span id="trim_mapping_li_20">trim_mapping_li_20.py</span></B><br/><!--代码位置 -->
		                     	<pre class="prettyprint linenums">
import os
import glob
import argparse
import subprocess
import multiprocessing


"""
this script is designed to trim and align RRBS/WGBS fastq files
with trim_galore and bismark
__author__ = likai
"""


def get_files(path_fastq):
    """
    directory to fastq files
    :param path_fastq:
    :return:
    """
    print(f"Start reading from {path_fastq}")
    fastq_files = glob.glob(f"{path_fastq}/*.fastq")
    return fastq_files


def trim_map_file(fastq, clean_dir, log_dir, bismark_index, bowtie2_path, mapping_dir, mapping_result_dir):
    """
    trim fastq files with trim_galore, and align them with bismark
    :param fastq: directory to fastq files
    :param clean_dir: directory to put trimmed fastq files
    :param log_dir: directory to put log files
    :param bismark_index: directory to bismark index files
    :param bowtie2_path: directory to bowtie2 executable
    :param mapping_dir: directory to put aligned bam files
    :param mapping_result_dir: directory to put bismark log files
    :return: none
    """
    prefix = os.path.basename(fastq).split(".")[0]
    this_sample_clean_out = f"{clean_dir}/{prefix}"
    this_sample_mapping_out = f"{mapping_dir}/{prefix}"
    this_sample_trim_log = f"{log_dir}/{prefix}.trim.log"
    this_sample_bismark_log = f"{mapping_result_dir}/{prefix}/{prefix}.bismark.txt"
    this_sample_error_log = f"{log_dir}/{prefix}.trim.error.log"
    if not os.path.exists(f"{mapping_result_dir}/{prefix}"):
        os.makedirs(f"{mapping_result_dir}/{prefix}")
    if not os.path.exists(this_sample_clean_out):
        os.makedirs(this_sample_clean_out)
    if not os.path.exists(this_sample_mapping_out):
        os.makedirs(this_sample_mapping_out)
    if fastq.endswith(".sra_1.fastq"):
        fastq_file_1 = f"{os.path.dirname(fastq)}/{prefix}.sra_1.fastq"
        fastq_file_2 = f"{os.path.dirname(fastq)}/{prefix}.sra_2.fastq"
        cmd_trim = f"trim_galore --rrbs --illumina --phred33 " \
                   f"--paired {fastq_file_1} {fastq_file_2} -o {this_sample_clean_out} > {this_sample_trim_log} 2>&1"
        clean_file_1 = f"{this_sample_clean_out}/{prefix}.sra_1_val_1.fq"
        clean_file_2 = f"{this_sample_clean_out}/{prefix}.sra_2_val_2.fq"
        cmd_bismark = f"bismark {bismark_index} --path_to_bowtie {bowtie2_path} --bowtie2 " \
                      f"-1 {clean_file_1} -2 {clean_file_2} -o {this_sample_mapping_out} " \
                      f"--temp_dir {log_dir} > {this_sample_bismark_log} 2>&1"
        bam_file = f"{this_sample_mapping_out}/{prefix}.sra_1_val_1_bismark_bt2_pe.bam"
    else:
        cmd_trim = f"trim_galore --rrbs --illumina --phred33 " \
                   f"{fastq} -o {this_sample_clean_out} > {this_sample_trim_log} 2>&1"
        clean_file = f"{this_sample_clean_out}/{prefix}.sra_trimmed.fq"
        cmd_bismark = f"bismark {bismark_index} --path_to_bowtie {bowtie2_path} " \
                      f"--bowtie2 {clean_file} -o {this_sample_mapping_out} " \
                      f"--temp_dir {log_dir} > {this_sample_bismark_log} 2>&1"
        bam_file = f"{this_sample_mapping_out}/{prefix}.sra_trimmed_bismark_bt2.bam"
    bam_qc_dir = f"{mapping_result_dir}/bam_qc/{prefix}"
    if not os.path.exists(bam_qc_dir):
        os.makedirs(bam_qc_dir)
    cmd_bam_qc = f"fastqc -o {bam_qc_dir} -q {bam_file}"
    try:
        print(f"  Trimming: {prefix}")
        subprocess.check_output(cmd_trim, shell=True)
        print(f"  Aligning: {prefix}")
        subprocess.check_output(cmd_bismark, shell=True)
        print(f"  Bam qc: {prefix}")
        subprocess.check_output(cmd_bam_qc, shell=True)
    except Exception as e:
        print(e)
        with open(this_sample_error_log, "w+") as out:
            out.write(str(e) + "\n")


def multi_process_run(fastq_files, num_process, clean_dir, log_dir, bismark_index, bowtie2_path, mapping_dir, mapping_result_dir):
    """
    to run faster, we create a multiprocessing pool here
    :param fastq_files: same as above, same below
    :param num_process:
    :param clean_dir:
    :param log_dir:
    :param bismark_index:
    :param bowtie2_path:
    :param mapping_dir:
    :param mapping_result_dir:
    :return:
    """
    pool = multiprocessing.Pool(processes=int(num_process))
    for fastq in fastq_files:
        if not fastq.endswith(".sra_2.fastq"):
            pool.apply_async(trim_map_file, (fastq, clean_dir, log_dir, bismark_index, bowtie2_path, mapping_dir, mapping_result_dir,))
    pool.close()
    pool.join()


def get_multiqc_file(mapping_dir, mapping_result_dir):
    """
    collect qc file for multiQC
    :param mapping_dir:
    :param mapping_result_dir:
    :return:
    """
    print(f"Start reading from {mapping_dir}")
    bismark_report_files = glob.glob(f"{mapping_dir}/*/*_report.txt")
    bismark_report_multiqc_path = f"{mapping_result_dir}/mapping_multiqc_sample.txt"
    print(f"Start reading from {mapping_result_dir}/bam_qc")
    bam_qc_files = glob.glob(f"{mapping_result_dir}/bam_qc/*/*.zip")
    bam_qc_multiqc_path = f"{mapping_result_dir}/bam_qc_multiqc_sample.txt"
    with open(bismark_report_multiqc_path, "w+") as out:
        out.write("\n".join(bismark_report_files) + "\n")
    with open(bam_qc_multiqc_path, "w+") as out:
        out.write("\n".join(bam_qc_files) + "\n")
    return bismark_report_multiqc_path, bam_qc_multiqc_path


def multiqc(bismark_report_multiqc_path, bam_qc_multiqc_path, n, mapping_result_dir):
    """
    run multiQC
    :param bismark_report_multiqc_path:
    :param bam_qc_multiqc_path:
    :param n:
    :param mapping_result_dir:
    :return:
    """
    print("Start multiqc")
    n_bismark = f"{n}_bismark"
    n_bam_qc = f"{n}_bam_qc"
    cmd_bismark = f"multiqc -n {n_bismark} -s -o {mapping_result_dir} -l {bismark_report_multiqc_path}"
    cmd_bam_qc = f"multiqc -n {n_bam_qc} -s -o {mapping_result_dir} -l {bam_qc_multiqc_path}"
    try:
        subprocess.check_output(cmd_bismark, shell=True)
        subprocess.check_output(cmd_bam_qc, shell=True)
    except Exception as e:
        print(e)


def main():
    parser = argparse.ArgumentParser(description="A program to map raw fastq and remove duplicates")
    parser.add_argument("-i", action="store", dest="path_fastq")
    parser.add_argument("-p", action="store", dest="num_process", help="processes to use in the program", default=10)
    parser.add_argument("-m", action="store", dest="mapping_result_dir", help="path to put mapping result in")
    parser.add_argument("-l", action="store", dest="log_dir", help="path to put log files")
    parser.add_argument("-b", action="store", dest="mapping_dir", help="path to put bam files")
    parser.add_argument("-n", action="store", dest="n", help="project name to use in MultiQC output[bioproject.13]")
    parser.add_argument("-t", action="store", dest="bowtie2_path", help="path to bowtie2 executable")
    parser.add_argument("-d", action="store", dest="path_index", help="path to bismark/bowtie2 indexes")
    parser.add_argument("-c", action="store", dest="clean_dir", help="path to put trimmed fastq files")
    results = parser.parse_args()
    if not all([results.path_fastq, results.num_process, results.mapping_result_dir, results.log_dir,
                results.n, results.mapping_dir, results.bowtie2_path, results.path_index, results.clean_dir]):
        print("too few arguments, type -h for more information!")
        exit(-1)
    path_fastq = results.path_fastq if not results.path_fastq.endswith("/") else results.path_fastq.rstrip("/")
    num_process = results.num_process
    mapping_result_dir = results.mapping_result_dir if not results.mapping_result_dir.endswith(
        "/") else results.mapping_result_dir.rstrip("/")
    log_dir = results.log_dir if not results.log_dir.endswith("/") else results.log_dir.rstrip("/")
    mapping_dir = results.mapping_dir if not results.mapping_dir.endswith("/") else results.mapping_dir.rstrip("/")
    n = results.n
    bowtie2_path = results.bowtie2_path if not results.bowtie2_path.endswith("/") else results.bowtie2_path.rstrip("/")
    path_index = results.path_index if not results.path_index.endswith("/") else results.path_index.rstrip("/")
    clean_dir = results.clean_dir if not results.clean_dir.endswith("/") else results.clean_dir.rstrip("/")
    for each_dir in [mapping_result_dir, log_dir, mapping_dir, clean_dir]:
        if not os.path.exists(each_dir):
            os.makedirs(each_dir)
    fastq_files = get_files(path_fastq)
    multi_process_run(fastq_files, num_process, clean_dir, log_dir,
                      path_index, bowtie2_path, mapping_dir, mapping_result_dir)
    bismark_report_multiqc_path, bam_qc_multiqc_path = get_multiqc_file(mapping_dir, mapping_result_dir)
    multiqc(bismark_report_multiqc_path, bam_qc_multiqc_path, n, mapping_result_dir)


if __name__ == '__main__':
    main()
		</pre>
		                     </p>
	                     </div>
                    </div> 
                  </div>
                </div> 
              </div>

            </div>
          </div>
        </div>
        <!-- /page content -->

        <!-- footer content -->
        <footer>
          <div class="pull-center">
            <strong>Copyright &copy; 2019 Big data Institute, Wenzhou Medical University.</strong> All rights reserved.<span><script type="text/javascript" src="//ra.revolvermaps.com/0/0/3.js?i=0e4hv6i8yu6&amp;b=0&amp;s=21&amp;m=0&amp;cl=ffffff&amp;co=010020&amp;cd=aa0000&amp;v0=60&amp;v1=60&amp;r=1" async="async"></script></span>
          </div>
          <div class="clearfix"></div>
        </footer>
        <!-- /footer content -->
      </div>
    </div>

    <!-- jQuery -->
    <script src="../js/jquery.min.js"></script>
    <!-- Bootstrap -->
    <script src="../js/bootstrap.min.js"></script>
    <!-- NProgress -->
    <script src="../js/nprogress.js"></script>
    <!-- Custom Theme Scripts -->
    <script src="../js/custom.min.js"></script> 
  </body>
</html>
Zuguang Gu