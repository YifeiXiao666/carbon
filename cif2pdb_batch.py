import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx
from biotite.structure.io.pdbx import CIFFile
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io import save_structure
import os
import gzip
import warnings
import glob
import traceback
import json

from biotite.structure.io.pdbx.convert import (
    _parse_intra_residue_bonds, 
    _parse_inter_residue_bonds, 
    _filter, 
    _get_block
)

PDB_ROOT = '/home/megagatlingpea/xiaoyifei/'

def parse_complete_cif(cif_path, dump_pdb_path=None, save_mode='all', sanity_check=True):

    sanity_failed = []
    if cif_path.endswith('.cif.gz'):
        cif_file = CIFFile.read(gzip.open(cif_path, "rt"))
    elif cif_path.endswith('.cif'):
        cif_file = CIFFile.read(cif_path)
    else:
        raise ValueError("Unsupported file format")


    structure_all_altlocs = pdbx.get_structure(
        cif_file, model=1, altloc='all', 
        extra_fields=['b_factor','occupancy','atom_id','charge'],
        include_bonds=False
    )

    block = _get_block(cif_file, None)
    atom_site_original = block.get("atom_site")
    struct_conn_original = block.get("struct_conn")
    chem_comp_bond_original = block.get("chem_comp_bond")

    if atom_site_original is None:
        raise ValueError("lack 'atom_site' in file, cannot proceed.")

    custom_bond_dict = None
    if chem_comp_bond_original:
        try:
            custom_bond_dict = _parse_intra_residue_bonds(chem_comp_bond_original)
        except KeyError:
            warnings.warn(
                "'_chem_comp_bond' lacks, use ccd.",
                UserWarning
            )

    altloc_ids = np.unique(structure_all_altlocs.altloc_id)
    conformation_ids = [alt_id for alt_id in altloc_ids if (alt_id != '.' and alt_id != ' ')]
    common_indices = np.where((structure_all_altlocs.altloc_id == '.') | (structure_all_altlocs.altloc_id == ' '))[0]

    all_bonds_list = []

    if not conformation_ids:
        bonds = struc.connect_via_residue_names(structure_all_altlocs, custom_bond_dict=custom_bond_dict)
        if struct_conn_original:
            inter_bonds = _parse_inter_residue_bonds(
                atom_site_original, struct_conn_original, atom_count=len(structure_all_altlocs)
            )
            bonds = bonds.merge(inter_bonds)
        structure_all_altlocs.bonds = bonds

    else:
        for conf_id in conformation_ids:
            conf_indices = np.where(structure_all_altlocs.altloc_id == conf_id)[0]
            temp_indices = np.concatenate([common_indices, conf_indices])
            temp_indices.sort() 
            
            temp_structure = structure_all_altlocs[temp_indices]

            temp_atom_site = _filter(atom_site_original, temp_indices)
            
            temp_bonds = struc.connect_via_residue_names(temp_structure, custom_bond_dict=custom_bond_dict)
            if struct_conn_original:
                inter_bonds = _parse_inter_residue_bonds(
                    temp_atom_site, struct_conn_original, atom_count=len(temp_structure)
                )
                temp_bonds = temp_bonds.merge(inter_bonds)

            original_indices_bonds = temp_bonds.as_array()
            if original_indices_bonds.shape[0] > 0:
                original_indices_bonds[:, 0] = temp_indices[original_indices_bonds[:, 0]]
                original_indices_bonds[:, 1] = temp_indices[original_indices_bonds[:, 1]]
                all_bonds_list.append(original_indices_bonds)

        if all_bonds_list:
            final_bond_array = np.vstack(all_bonds_list)
            final_bond_array[:, :2] = np.sort(final_bond_array[:, :2], axis=1)
            _, unique_indices = np.unique(final_bond_array, axis=0, return_index=True)
            unique_bonds = final_bond_array[np.sort(unique_indices)]
            
            final_bond_list = struc.BondList(structure_all_altlocs.array_length(), bonds=unique_bonds)
            structure_all_altlocs.bonds = final_bond_list

    if sanity_check:
        
        fname = os.path.basename(cif_path)
        pdb_id = fname.replace('.cif.gz', '').replace('.cif', '')
        pdb_file = f'/mnt/data/wwpdb/data/structures/all/pdb/pdb{pdb_id}.ent.gz'
        pdb_file = PDBFile.read(gzip.open(pdb_file, "rt"))

        cif_struc_first = pdbx.get_structure(
            cif_file, model=1, altloc='first', 
            extra_fields=['b_factor','occupancy','atom_id','charge'],
            include_bonds=True
        )

        pdb_struc_first = pdb.get_structure(
            pdb_file, model=1, altloc='first', 
            extra_fields=['b_factor','occupancy','atom_id','charge'],
            include_bonds=True
        )

        if len(conformation_ids) > 0:
            flag = []
            for alt in conformation_ids:
                structure_all_altlocs_model_alt = structure_all_altlocs[(structure_all_altlocs.altloc_id == '.') | (structure_all_altlocs.altloc_id == ' ') | (structure_all_altlocs.altloc_id == alt)]
                try:
                    if ((structure_all_altlocs_model_alt.bonds._bonds - cif_struc_first.bonds._bonds).max() < 1e-6 ) and ((structure_all_altlocs_model_alt.coord - cif_struc_first.coord).max() < 1e-6) and (structure_all_altlocs_model_alt.bonds._bonds.shape == pdb_struc_first.bonds._bonds.shape):
                        flag.append(alt)
                    else:
                        continue
                except:
                    continue
                    
                
            
            if len(flag) == 0:
                raise Exception('Sanity check failed')

        else:
            # sanitizing all altlocs is the same with first altloc selection when no altlocs
            try:
                assert (structure_all_altlocs.bonds._bonds - cif_struc_first.bonds._bonds).max() < 1e-6
                assert (structure_all_altlocs.coord - cif_struc_first.coord).max() < 1e-6
                assert structure_all_altlocs.bonds._bonds.shape == pdb_struc_first.bonds._bonds.shape
            except AssertionError:
                sanity_failed.append(cif_path.split('/')[-1][0:4])
                print(f"Sanity check failed for {cif_path.split('/')[-1][0:4]}.")

    if dump_pdb_path is not None:
        if save_mode == 'first':
            save_structure(dump_pdb_path, cif_struc_first)
        elif save_mode == 'all':
            save_structure(dump_pdb_path, structure_all_altlocs)
        elif save_mode == 'separate':
            if len(conformation_ids) == 0:
                save_structure(dump_pdb_path, cif_struc_first)
            else:
                # for alt_id in conformation_ids:
                for alt_id in flag:
                    conf_structure = structure_all_altlocs[(structure_all_altlocs.altloc_id == '.') | (structure_all_altlocs.altloc_id == ' ') | (structure_all_altlocs.altloc_id == alt_id)]
                    conf_dump_path = dump_pdb_path.replace('.pdb', f'_altloc_{alt_id}.pdb')
                    save_structure(conf_dump_path, conf_structure)
        else:
            raise ValueError("Unsupported save_mode. Choose from 'first', 'all', or 'separate'.")

    return structure_all_altlocs, sanity_failed

def parse_complete_cif_batch(input_path, dump_root=None, save_mode='all', sanity_check=True, pattern='*.cif.gz'):
    """
    Batch-convert cif files to pdb using parse_complete_cif.

    Parameters:
    - input_path: directory path, single file path, or list of file paths
    - dump_root: directory to write pdb outputs (optional). If None, no files are written unless parse_complete_cif was called with a dump path.
    - save_mode: 'first'|'all'|'separate' passed to parse_complete_cif
    - sanity_check: whether to run sanity checks
    - pattern: glob pattern when input is a directory (default '*.cif.gz')
    Returns a summary dict with keys: processed, success, failed.
    """
    results = {"processed": 0, "success": [], "failed": []}
    # normalize input to file list
    if isinstance(input_path, (list, tuple)):
        files = list(input_path)
    else:
        if os.path.isfile(input_path):
            files = [input_path]
        else:
            files = sorted(glob.glob(os.path.join(input_path, pattern)))

    for cif in files:
        try:
            fname = os.path.basename(cif)
            base = fname.replace('.cif.gz','').replace('.cif','')
            dump_pdb = None
            if dump_root is not None:
                os.makedirs(dump_root, exist_ok=True)
                dump_pdb = os.path.join(dump_root, base + '.pdb')
            # call existing function
            structure, sanity_failed = parse_complete_cif(cif, dump_pdb_path=dump_pdb, save_mode=save_mode, sanity_check=sanity_check)
            results['processed'] += 1
            if sanity_failed:
                results['failed'].append((cif, sanity_failed))
            else:
                results['success'].append(cif)
        except Exception as e:
            results['processed'] += 1
            tb = traceback.format_exc()
            results['failed'].append((cif, str(e) + '\n' + tb))
    return results

# Example usage (uncomment and edit paths to run):
summary = parse_complete_cif_batch('/home/megagatlingpea/xiaoyifei/debug/', dump_root='/home/megagatlingpea/xiaoyifei/pdb_dumps', save_mode='separate', sanity_check=True)
print(json.dumps(summary, indent=2))