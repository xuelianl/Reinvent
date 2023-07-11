# /pubhome/xli02/opt/miniconda/envs/reinvent.v3.2/lib/python3.7/site-packages/reinvent_scoring/scoring/score_components/structural/shape_similarity.py

import os
import json
import shutil
import tempfile
import numpy as np
import pandas as pd
from typing import List

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_summary import ComponentSummary
# from reinvent_scoring.scoring.score_components import BaseScoreComponent
from reinvent_scoring.scoring.score_components.structural.base_structural_component import BaseStructuralComponent


class ShapeSimilarity(BaseStructuralComponent):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)
        self._configuration_path = self.parameters.specific_parameters[self.component_specific_parameters.SIMI_CONFIG_PATH]
        if isinstance(self._configuration_path, str):
            if os.path.isfile(self._configuration_path):
                with open(self._configuration_path) as file:
                    conf = file.read().replace("\r", "").replace("\n", "")
            conf = json.loads(conf)
        self._conf = conf

    # def calculate_score(self, molecules: List) -> ComponentSummary:
    #     score = self._calculate_qed(molecules)
    #     score_summary = ComponentSummary(total_score=score, parameters=self.parameters)
    #     return score_summary

    def _calculate_3D_simi_shafts(self, poses_path, scores_path, pref, template_idx, template_3d_dir):
        def apply_prefix_to_filename(path: str, output_prefix: str):
            """This method modifies the output path directory by concatenating an output_prefix at the front

            :param path:: The base output path pre-modification
            :type path: string
            :param output_prefix: prefix to be concatenated at the front to modify the output path
            :type output_prefix: string
            :return: string path with output_prefix concatenated at the front if path is not None.
                Otherwise, return the base path
            """
            if output_prefix is not None:
                return os.path.join(os.path.dirname(path), "".join([output_prefix, os.path.basename(path)]))
            else:
                return path

        # note that "text" is True (in contrast to the underlying "mkstemp")
        def gen_temp_file(suffix=None, prefix=None, dir=None, text=True) -> str:
            filehandler, path = tempfile.mkstemp(suffix=suffix, prefix=prefix, dir=dir, text=text)
            os.close(filehandler)
            return path
        
        poses_path = apply_prefix_to_filename(poses_path, pref)

        scores_path = apply_prefix_to_filename(scores_path, pref)
        if not os.path.exists(scores_path):
            return pd.DataFrame()
        score_df = pd.read_csv(scores_path)

        pose_mol2 = gen_temp_file(suffix=".mol2")
        # pose_mol2 = '/tmp/' + "".join([pref, os.path.basename(poses_path).split('.')[0]]).replace('"', '') + '.mol2'
        os.system(f'/pubhome/xli02/Downloads/ZBH/unicon/unicon -i {poses_path} -o {pose_mol2} -v 0')

        # inception_df = pd.read_csv('/pubhome/xli02/project/mpro/reinvent/inception_from_structures/unique_structures_filtered_N_504.csv', sep='\t') #
        inception_df = pd.read_csv(template_idx, sep='\t')
        # template_3d_dir = '/pubhome/xli02/project/mpro/reinvent/3-3D_shape_simi/inception_502_aligned_7VTH_lig_mol2' #
        incep_tmp_dir = tempfile.mkdtemp()
        for incep in inception_df['crystal_name']:
            incep_mol2 = f'{template_3d_dir}/incep_{incep}.mol2'
            os.system(f'mkdir {incep_tmp_dir}/shafts_{incep}; cd {incep_tmp_dir}/shafts_{incep}; Cynthia -q {incep_mol2} -t {pose_mol2} -scoreOnly &> shafts.log')
        hybrid_df = score_df[['name']].copy()
        # shape_df = score_df[['name']]
        # feature_df = score_df[['name']]
        for incep in inception_df['crystal_name']:
            res = pd.read_csv(f'{incep_tmp_dir}/shafts_{incep}/Result.list', sep='\t')
            hybrid_df = hybrid_df.merge(res[['Name', 'HybridScore']].copy().rename(columns={'HybridScore':incep, 'Name':'name'}), on='name', how='left')
            # shape_df = shape_df.merge(res[['Name', 'ShapeScore']].copy().rename(columns={'ShapeScore':incep, 'Name':'name'}), on='name', how='left')
            # feature_df = feature_df.merge(res[['Name', 'FeatureScore']].copy().rename(columns={'FeatureScore':incep, 'Name':'name'}), on='name', how='left')
            shutil.rmtree(f'{incep_tmp_dir}/shafts_{incep}')
        shutil.rmtree(incep_tmp_dir)
        if os.path.exists(pose_mol2):
            os.remove(pose_mol2)
        hybrid_df_set_idx = hybrid_df.set_index('name')
        hybrid_max_df = hybrid_df_set_idx.idxmax(axis="columns").reset_index().rename(columns={0: 'max_hybrid_cry'}).merge(hybrid_df_set_idx.max(axis=1).reset_index().rename(columns={0:'max_hybrid'}), on = 'name')
        # simi_dir = os.path.join(os.path.dirname(os.path.dirname(poses_path)), 'shafts_simi')
        # if not os.path.exists(simi_dir):
        #     os.makedirs(simi_dir)
        # hybrid_max_df.to_csv(simi_dir + '/' + pref + 'shafts_hybrid.csv', sep='\t', index=False)
        hybrid_max_df['idx'] = hybrid_max_df['name'].str.split(':').str[0].astype(int)
        best_simi_df = hybrid_max_df.groupby('idx')['max_hybrid'].agg('max').reset_index()
        return best_simi_df


    # def calculate_3D_simi_rdkit(self, docking_run, docker, output_prefix):
    #     from collections import defaultdict
    #     from rdkit import Chem
    #     from rdkit.Chem import AllChem, rdShapeHelpers
    #     poses_path = docking_run[_DE.OUTPUT][_DE.OUTPUT_POSES][_DE.OUTPUT_POSES_PATH]
    #     poses_path = docker.apply_prefix_to_filename(poses_path, output_prefix)
    #     pref = output_prefix.replace('"', '')
    #     logger.log(f'{output_prefix}, {pref}', _LE.INFO) #
    #     logger.log(poses_path, _LE.INFO) #
        
    #     with open(poses_path, 'r') as f:
    #         lines = f.readlines()
        
    #     mol_start = True
    #     docked_idx2sdf_lines = defaultdict(list)
    #     for line in lines:
    #         if mol_start:
    #             mol_idx = line.rstrip('\n').replace(':', '_')
    #         if '$$$$' in line:
    #             docked_idx2sdf_lines[mol_idx].append('$$$$\n')
    #             mol_start = True
    #         else:
    #             mol_start = False
    #             docked_idx2sdf_lines[mol_idx].append(line)
        
    #     inception_df = pd.read_csv('/pubhome/xli02/project/mpro/reinvent/inception_from_structures/unique_structures_filtered_N_504.csv', sep='\t') #
    #     aligned_dir = '/pubhome/xli02/project/mpro/reinvent/aligned_stru/7VTH/all_504_inception_20230302_for_3D_simi/aligned_structures' #
    #     cry_ligs_to_pdb_lines = defaultdict(list)
    #     for row in inception_df.itertuples():
    #         cry_pdb = f'{aligned_dir}/{row.crystal_name}_aligned.pdb'
            
    #         with open(cry_pdb, 'r') as f2:
    #             cry_pdb_lines = f2.readlines()
    #         for line in cry_pdb_lines:
    #             if line[:6] == 'HETATM':
    #                 if len(row.lig_ident) == 9 and line[17:26]==row.lig_ident:
    #                     cry_ligs_to_pdb_lines[row.crystal_name].append(line)
    #                 elif len(row.lig_ident) == 10 and line[16:26]==row.lig_ident:
    #                     cry_ligs_to_pdb_lines[row.crystal_name].append(line)
        
    #     best_list = []
    #     failed_pose = []
    #     failed_cry = []
    #     for mol_idx in docked_idx2sdf_lines.keys():
    #         # pose_sdf = f'/tmp/{output_prefix}m{mol_idx}_docked_poses.sdf'
    #         pose_sdf = f'/tmp/{pref}m{mol_idx}_docked_poses.sdf'
    #         with open(pose_sdf, 'w') as f3:
    #             for line in docked_idx2sdf_lines[mol_idx]:
    #                 f3.write(line)
    #         if Chem.SDMolSupplier(pose_sdf)[0] is None:
    #             failed_pose.append(mol_idx)
    #             best_list.append([mol_idx.replace('_', ':'), 1])
    #             if os.path.exists(pose_sdf):
    #                 os.remove(pose_sdf)
    #             continue
    #         else:
    #             pose_sdf_mol = Chem.SDMolSupplier(pose_sdf)[0]
    #         if os.path.exists(pose_sdf):
    #             os.remove(pose_sdf)
    #         std = []
    #         for cry in cry_ligs_to_pdb_lines.keys():
    #             cry_mol = Chem.MolFromPDBBlock(''.join(cry_ligs_to_pdb_lines[cry]))
    #             if cry_mol is not None and cry_mol.GetNumAtoms() != 0:
    #                 std.append(rdShapeHelpers.ShapeTanimotoDist(pose_sdf_mol, cry_mol))
    #             else:
    #                 if cry not in failed_cry:
    #                     failed_cry.append(cry)
    #         best_list.append([mol_idx.replace('_', ':'), min(std)])
    #     best_simi_df = pd.DataFrame(best_list, columns=['docked_pose_idx', 'min_std'])
    #     shape_similarity_dir = f'{os.path.dirname(os.path.dirname(poses_path))}/shape_similarity'
    #     if not os.path.exists(shape_similarity_dir):
    #         os.makedirs(shape_similarity_dir)
    #     best_simi_df.to_csv(f'{shape_similarity_dir}/{pref}shape_simi.csv', sep='\t', index=False)
    #     logger.log(f'rdkit_read_failed_cry when calculating 3D similarity: {failed_cry}.', _LE.INFO) #
    #     logger.log(f'rdkit_read_failed_sdf_idx when calculating 3D similarity: {failed_pose}', _LE.INFO) #
    #     return best_simi_df

    def _calculate_score(self, query_mols, step) -> np.array:
        output_prefix = self._get_step_string(step)
        pref = output_prefix.replace('"', '')
        # best_simi_df = pd.read_csv(f'{shape_similarity_dir}/{pref}shape_simi.csv', sep='\t')
        poses_path = self._conf['similarity_3d']['poses_path']
        scores_path = self._conf['similarity_3d']['scores_path']
        template_idx = self._conf['similarity_3d']['template_idx']
        template_3d_dir = self._conf['similarity_3d']['template_3d_dir']
        if self._conf['similarity_3d']['simi_method'] == 'shafts_hybrid':
            best_simi_df = self._calculate_3D_simi_shafts(poses_path, scores_path, pref, template_idx, template_3d_dir)
        if len(best_simi_df) == 0:
            return np.zeros(len(query_mols)), np.zeros(len(query_mols))
        # best_simi_df['idx'] = best_simi_df['name'].str.split(':').str[0].astype(int)
        query_mol_idx_df = pd.DataFrame(range(len(query_mols)), columns=['idx'])
        best_simi_df_merged = query_mol_idx_df.merge(best_simi_df, on='idx', how='left')
        scores = []
        for score in best_simi_df_merged['max_hybrid']:
            if np.isnan(score):
                score = 0
            scores.append(score)

        # template_3d_dir = self._conf['similarity_3d']['template_3d_dir']
        transform_params = self.parameters.specific_parameters.get(
            self.component_specific_parameters.TRANSFORMATION, {}
        )
        transformed_scores = self._transformation_function(scores, transform_params)

        return np.array(transformed_scores), np.array(scores)
        # return np.array(1-best_simi_df['min_std'], dtype=np.float32)

    def _parse_result(self, result):
        pass

    def _create_command(self, smiles: List[str], step):
        pass
