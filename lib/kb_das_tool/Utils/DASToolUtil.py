import errno
import json
import os
import subprocess
import sys
import time
import uuid
import zipfile

from Bio import SeqIO

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils
from installed_clients.ReadsUtilsClient import ReadsUtils


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class DASToolUtil:
    DASTOOL_THREADS=2

    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.shock_url = config['shock-url']
        self.ws_url = config['workspace-url']
        self.dfu = DataFileUtil(self.callback_url)
        self.ru = ReadsUtils(self.callback_url)
        self.au = AssemblyUtil(self.callback_url)
        self.mgu = MetagenomeUtils(self.callback_url)



    def validate_run_das_tool_params(self, params):
        """
        validate_run_concoct_params:
                validates params passed to run_concoct method

        """
        log('Start validating run_kb_das_tool params')

        # check for required parameters
        for p in ['binned_contig_name', 'workspace_name']:
            if p not in params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

    def mkdir_p(self, path):
        """
        mkdir_p: make directory for given path
        """
        if not path:
            return
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def run_command(self, command):
        """
        run_command: run command and print result
        """
        os.chdir(self.scratch)
        log('Start executing command:\n{}'.format(command))
        log('Command is running from:\n{}'.format(self.scratch))
        pipe = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        output,stderr = pipe.communicate()
        exitCode = pipe.returncode

        if (exitCode == 0):
            log('Executed command:\n{}\n'.format(command) +
                'Exit Code: {}\n'.format(exitCode))
        else:
            error_msg = 'Error running command:\n{}\n'.format(command)
            error_msg += 'Exit Code: {}\nOutput:\n{}\nStderr:\n{}'.format(exitCode, output, stderr)
            raise ValueError(error_msg)
            sys.exit(1)
        return (output,stderr)



    def get_contig_file(self, assembly_ref):
        """
        get_contig_file: get contif file from GenomeAssembly object
        """

        contig_file = self.au.get_assembly_as_fasta({'ref': assembly_ref}).get('path')

        sys.stdout.flush()
        contig_file = self.dfu.unpack_file({'file_path': contig_file})['file_path']

        return contig_file


    def retrieve_and_clean_assembly(self, task_params):
        if os.path.exists(task_params['contig_file_path']):
            assembly =  task_params['contig_file_path']
            print("FOUND ASSEMBLY ON LOCAL SCRATCH")
        else:
            # we are on njsw so lets copy it over to scratch
            assembly = self.get_contig_file(task_params['assembly_ref'])

        # remove spaces from fasta headers because that breaks bedtools
        assembly_clean = os.path.abspath(assembly).split('.fa')[0] + "_clean.fa"

        command = '/bin/bash reformat.sh in={} out={} addunderscore'.format(assembly,assembly_clean)

        log('running reformat command: {}'.format(command))
        out,err = self.run_command(command)

        return assembly_clean


    def generate_das_tool_command(self, params):
        """
        generate_command: generate concoct params
        """
        # print("params" + str(params))
        # params['contig_file_path'] = self.retrieve_and_clean_assembly(params)
        # print("params" + str(params))

        binned_contig_files = self.generate_input_binned_contig_list_command[0]
        binned_contig_names = self.generate_input_binned_contig_list_command[1]


        print("\n\nRunning generate_das_tool_command")

        command = 'DAS_Tool '

        command += '-i {} '.format(binned_contig_files)
        command += '-l {} '.format(binned_contig_names)
        command += '-c {} '.format(params.get('contig_file_path'))
        command += '-o das_tool_output_dir '
        command += '--search_engine {} '.format(params.get('search_engine'))
        command += '--score_threshold {} '.format(params.get('score_threshold'))
        command += '--duplicate_penalty {} '.format(params.get('duplicate_penalty'))
        command += '--megabin_penalty {} '.format(params.get('megabin_penalty'))
        command += '--write_bin_evals {} '.format(params.get('write_bin_evals'))
        command += '--create_plots {} '.format(params.get('create_plots'))
        command += '--write_bins {} '.format(params.get('write_bins'))
        command += '--write_unbinned {} '.format(params.get('write_unbinned'))
        command += '-t {}'.format(self.DASTOOL_THREADS)

        log('Generated das_tool command: {}'.format(command))

        return command

    def generate_output_file_list(self, result_directory):
        """
        generate_output_file_list: zip result files and generate file_links for report
        """
        log('Start packing result files')
        output_files = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self.mkdir_p(output_directory)
        result_file = os.path.join(output_directory, 'das_tool_result.zip')
        report_file = None

        with zipfile.ZipFile(result_file, 'w',
                             zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zip_file:

            # grab all files we want to zip
            for dirname, subdirs, files in os.walk(result_directory):
                if (dirname.endswith(self.CONCOCT_BIN_DIR)):
                    baseDir=os.path.basename(dirname)
                    for file in files:
                        if (file.endswith('.summary')):
                            continue
                        full=os.path.join(dirname, file)
                        zip_file.write(full,os.path.join(baseDir,file))
                for file in files:
                    if (file.endswith('.depth.txt')):
                        zip_file.write(os.path.join(dirname, file),file)


        output_files.append({'path': result_file,
                             'name': os.path.basename(result_file),
                             'label': os.path.basename(result_file),
                             'description': 'Files generated by kb_das_tool App'})

        return output_files

    def generate_html_report(self, result_directory, assembly_ref, binned_contig_obj_ref, header):
        """
        generate_html_report: generate html summary report
        """

        log('Start generating html report')
        html_report = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self.mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'report.html')

        # get summary data from existing assembly object and bins_objects
        Summary_Table_Content = ''
        Overview_Content = ''
        (binned_contig_count, input_contig_count,
         total_bins_count) = self.generate_overview_info(assembly_ref,
                                                          binned_contig_obj_ref,
                                                          result_directory)

        Overview_Content += '<p>Binned contigs: {}</p>'.format(binned_contig_count)
        Overview_Content += '<p>Input contigs: {}</p>'.format(input_contig_count)
        Overview_Content += '<p>Number of bins: {}</p>'.format(total_bins_count)

        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'report_template.html'),
                      'r') as report_template_file:
                report_template = report_template_file.read()
                report_template = report_template.replace('<p>Overview_Content</p>',
                                                          Overview_Content)
                report_template = report_template.replace('Summary_Table_Content',
                                                          Summary_Table_Content)
                result_file.write(report_template)

        html_report.append({'path': result_file_path,
                            'name': os.path.basename(result_file_path),
                            'label': os.path.basename(result_file_path),
                            'description': 'HTML summary report for kb_concoct App'})
        return html_report


    def generate_overview_info(self, assembly_ref, binned_contig_obj_ref, result_directory):
        """
        _generate_overview_info: generate overview information from assembly and binnedcontig
        """

        # get assembly and binned_contig objects that already have some data populated in them
        assembly = self.dfu.get_objects({'object_refs': [assembly_ref]})['data'][0]
        binned_contig = self.dfu.get_objects({'object_refs': [binned_contig_obj_ref]})['data'][0]

        input_contig_count = assembly.get('data').get('num_contigs')
        bins_directory = os.path.join(self.scratch, result_directory, self.CONCOCT_BIN_DIR)
        binned_contig_count = 0
        total_bins_count = 0
        total_bins = binned_contig.get('data').get('bins')
        total_bins_count = len(total_bins)
        for bin in total_bins:
            binned_contig_count += len(bin.get('contigs'))

        return (binned_contig_count, input_contig_count, total_bins_count)


    def generate_report(self, binned_contig_obj_ref, params):
        """
        generate_report: generate summary report

        """
        log('Generating report')

        output_files = self.generate_output_file_list(params['result_directory'])

        output_html_files = self.generate_html_report(params['result_directory'],
                                                       params['assembly_ref'],
                                                       binned_contig_obj_ref)

        report_params = {
              'message': '',
              'workspace_name': params.get('workspace_name'),
              'file_links': output_files,
              'html_links': output_html_files,
              'direct_html_link_index': 0,
              'html_window_height': 266,
              'report_object_name': 'kb_das_tool_report_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def run_das_tool(self, params):
        """
        run_das_tool: DAS_Tool app

        required params:
            assembly_ref: Metagenome assembly object reference
            binned_contig_name: list of BinnedContig objects
            workspace_name: the name of the workspace it gets saved to.

        optional params:
            min_contig_length: minimum contig length; default 1000

        ref: https://github.com/cmks/DAS_Tool
        """
        log('--->\nrunning ConcoctUtil.run_concoct\n' +
            'params:\n{}'.format(json.dumps(params, indent=1)))

        self.validate_run_das_tool_params(params)

        print("\n\nFinished running validate_run_das_tool_params")
        #
        # contig_file = self.get_contig_file(params.get('assembly_ref'))
        # params['contig_file_path'] = contig_file

        #params['binned_contig_list_file'] = binned_contig_list_file

        #self.mgu.binned_contigs_to_file(params)

        result_directory = os.path.join(self.scratch, str("dastool_output_dir"))

        params['result_directory'] = result_directory

        self.mkdir_p(result_directory)

        cwd = os.getcwd()
        log('changing working dir to {}'.format(result_directory))
        os.chdir(result_directory)


        #run concoct
        # command = self.generate_das_tool_command(params)
        #
        # self.run_command(command)

        os.chdir(cwd)
        log('changing working dir to {}'.format(cwd))

        log('Saved result files to: {}'.format(result_directory))
        log('Generated files:\n{}'.format('\n'.join(os.listdir(result_directory))))
        #
        # generate_binned_contig_param = {
        #     'file_directory': result_directory,
        #     'assembly_ref': params.get('assembly_ref'),
        #     'binned_contig_name': params.get('binned_contig_name'),
        #     'workspace_name': params.get('workspace_name')
        # }
        #
        # binned_contig_obj_ref = self.mgu.file_to_binned_contigs(
        #                             generate_binned_contig_param).get('binned_contig_obj_ref')
        #
        # reportVal = self.generate_report(binned_contig_obj_ref, result_directory, params)

        binned_contig_obj_ref = 'test'
        returnVal = {
            'result_directory': result_directory,
            'binned_contig_obj_ref': binned_contig_obj_ref
        }

        #returnVal.update(reportVal)

        return returnVal
