# -*- coding: utf-8 -*-
import os
import time
import unittest
import shutil
from configparser import ConfigParser

from kb_das_tool.kb_das_toolImpl import kb_das_tool
from kb_das_tool.kb_das_toolServer import MethodContext
from kb_das_tool.authclient import KBaseAuth as _KBaseAuth

from kb_das_tool.Utils.DASToolUtil import DASToolUtil

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.ReadsUtilsClient import ReadsUtils
from installed_clients.WorkspaceClient import Workspace
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils


class kb_das_toolTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_das_tool'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(cls.token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': cls.token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_das_tool',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = kb_das_tool(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_kb_das_tool_" + str(suffix)
        # ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

        cls.ws_info = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa
        cls.dfu = DataFileUtil(os.environ['SDK_CALLBACK_URL'], token=cls.token)
        cls.ru = ReadsUtils(os.environ['SDK_CALLBACK_URL'], token=cls.token)
        cls.au = AssemblyUtil(os.environ['SDK_CALLBACK_URL'], token=cls.token)
        cls.mgu = MetagenomeUtils(os.environ['SDK_CALLBACK_URL'], token=cls.token)
        cls.das_tool_runner = DASToolUtil(cls.cfg)
        cls.prepare_data()

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    @classmethod
    def prepare_data(cls):
        """
        Lets put everything on workspace
        """
        #

    #     # READS 2
    #     # building Interleaved library
    #     pe2_reads_filename = 'lib2.oldstyle.fastq'
    #     pe2_reads_path = os.path.join(cls.scratch, pe2_reads_filename)
    #
    #     # gets put on scratch. "work/tmp" is scratch
    #     shutil.copy(os.path.join("data", pe2_reads_filename), pe2_reads_path)
    #
    #     int2_reads_params = {
    #         'fwd_file': pe2_reads_path,
    #         'sequencing_tech': 'Unknown',
    #         'wsname': cls.ws_info[1],
    #         'name': 'MyInterleavedLibrary2',
    #         'interleaved': 'true'
    #     }
    #
    #     #from scratch upload to workspace
    #     cls.int2_oldstyle_reads_ref = cls.ru.upload_reads(int2_reads_params)['obj_ref']
    #     #
        # building Assembly
        #
        assembly_filename = 'small_arctic_assembly.fa'
        cls.assembly_filename_path = os.path.join(cls.scratch, assembly_filename)
        shutil.copy(os.path.join("data", assembly_filename), cls.assembly_filename_path)

        # from scratch upload to workspace
        assembly_params = {
            'file': {'path': cls.assembly_filename_path},
            'workspace_name': cls.ws_info[1],
            'assembly_name': 'MyAssembly'
        }

        # puts assembly object onto shock
        cls.assembly_ref = cls.au.save_assembly_from_fasta(assembly_params)

        print('\nDone uploading assembly')

        # Genome Bin set 1
        genome_bin_folder_name = 'bins_concoct'
        print("\n\n\n\ngenome_bin_folder_name: {}".format(genome_bin_folder_name))
        genome_bin_folder_path = os.path.join(cls.scratch, genome_bin_folder_name)
        print("\n\n\n\ngenome_bin_folder_path: {}".format(genome_bin_folder_path))
        os.listdir(genome_bin_folder_path)

        task_params = []
        for dirname, subdirs, files in os.walk(genome_bin_folder_path):
            for file in files:
                print("file: {}".format(file))
                if file.endswith('.fasta'):
                    task_params['result_directory'] = os.path.join(self.scratch, str("dastool_output_dir"))
                    task_params['bin_result_directory'] = genome_bin_folder_name
                    das_tool_runner.make_binned_contig_summary_file_for_binning_apps(task_params)
        # gets put on scratch. "work/tmp" is scratch
        shutil.copytree(os.path.join("data", genome_bin_folder_name), genome_bin_folder_path)

        binned_contig_object_params = {
            'file_directory': genome_bin_folder_path,
            'assembly_ref': cls.assembly_ref,
            'binned_contig_name': 'concoct.test_data.BinnedContig',
            'workspace_name': cls.ws_info[1],
        }

        #from scratch upload to workspace

        cls.concoct_genome_bin_ref = cls.mgu.file_to_binned_contigs(binned_contig_object_params)
        print("cls.concoct_genome_bin_ref***")
        print(str(cls.concoct_genome_bin_ref))




        # Genome Bin set 2
        genome_bin_folder_name = 'bins_metabat'
        genome_bin_folder_path = os.path.join(cls.scratch, genome_bin_folder_name)
        print("genome_bin_folder_path: {}".format(genome_bin_folder_path))

        # gets put on scratch. "work/tmp" is scratch
        shutil.copytree(os.path.join("data", genome_bin_folder_name), genome_bin_folder_path)

        binned_contig_object_params = {
            'file_directory': genome_bin_folder_path,
            'assembly_ref': cls.assembly_ref,
            'binned_contig_name': 'metabat.test_data.BinnedContig',
            'workspace_name': cls.ws_info[1],
        }

        #from scratch upload to workspace

        cls.metabat_genome_bin_ref = cls.mgu.file_to_binned_contigs(binned_contig_object_params)
        print("cls.metabat_genome_bin_ref***")
        print(str(cls.metabat_genome_bin_ref))


        # Genome Bin set 3
        genome_bin_folder_name = 'bins_maxbin'
        genome_bin_folder_path = os.path.join(cls.scratch, genome_bin_folder_name)
        print("genome_bin_folder_path: {}".format(genome_bin_folder_path))

        # gets put on scratch. "work/tmp" is scratch
        shutil.copytree(os.path.join("data", genome_bin_folder_name), genome_bin_folder_path)

        binned_contig_object_params = {
            'file_directory': genome_bin_folder_path,
            'assembly_ref': cls.assembly_ref,
            'binned_contig_name': 'maxbin.test_data.BinnedContig',
            'workspace_name': cls.ws_info[1],
        }

        #from scratch upload to workspace

        cls.maxbin_genome_bin_ref = cls.mgu.file_to_binned_contigs(binned_contig_object_params)
        print("cls.maxbin_genome_bin_ref***")
        print(str(cls.maxbin_genome_bin_ref))


    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        return self.ws_info[1]

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def test_your_method(self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods
        ret = self.serviceImpl.run_kb_das_tool(self.ctx, {'assembly_ref': self.assembly_ref, 'workspace_name': self.wsName,
                    'binned_contig_names': [self.concoct_genome_bin_ref, self.metabat_genome_bin_ref, self.maxbin_genome_bin_ref] })
