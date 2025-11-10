import argparse
from pathlib import Path
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx
import os
import math
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

def convert_single_cif_to_pdb(args):
    """转换单个CIF文件到PDB"""
    cif_path, output_dir = args
    pdb_path = output_dir / (cif_path.stem + ".pdb")
    
    try:
        cif_file = pdbx.PDBxFile.read(str(cif_path))
        structure = pdbx.get_structure(
            cif_file, model=1, altloc='first', 
            extra_fields=['b_factor','occupancy','atom_id','charge'],
            include_bonds=True
        )

        pdb_file = pdb.PDBFile()
        pdb_file.set_structure(structure)
        pdb_file.write(str(pdb_path))

        return (True, cif_path, pdb_path)
    except Exception as e:
        return (False, cif_path, str(e))

def convert_cif_to_pdb(input_dir, output_dir, num_processes=None):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 确定进程数量
    if num_processes is None:
        # 默认使用CPU数量的80%
        num_processes = 128 # we fixed 128 for docker cannot well recognize cpu by code above.
        print(f"使用 {num_processes} 个进程 (CPU数量的80%)")

    # boltzdesign is different from Odesign sample
    if any(input_dir.rglob("*wounresol.cif")):
        # results come from Odesign
        general = "*wounresol.cif"
    else:
        # results come from boltzdesign
        general = "*.cif"

    # 收集所有CIF文件路径
    cif_files = list(input_dir.rglob(general))
    total_files = len(cif_files)
    print(f"找到 {total_files} 个CIF文件待转换")
    
    total_files = len(cif_files)
    if total_files == 0:
        print("⚠️ 警告: 没有找到以指定前缀开头的CIF文件")
        print(f"查找的文件名前缀: {', '.join(target_prefixes)}")
    else:
        print(f"找到 {total_files} 个符合条件的CIF文件待转换")

    # 准备参数列表
    args_list = [(cif_path, output_dir) for cif_path in cif_files]

    # 使用进程池并行处理
    with Pool(num_processes) as pool:
        results = list(tqdm(
            pool.imap(convert_single_cif_to_pdb, args_list),
            total=total_files,
            desc="转换进度"
        ))
    
    # 统计结果
    success = sum(1 for r in results if r[0])
    failed = total_files - success
    
    # 打印失败的转换（如果有）
    if failed > 0:
        print("\n失败的转换:")
        for r in results:
            if not r[0]:
                print(f"❌ {r[1]} → 错误: {r[2]}")
    
    print(f"\n转换完成: {success}/{total_files} 成功, {failed} 失败")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="递归将 CIF 转换为 PDB 并保留目录结构 (多进程版)")
    parser.add_argument("--input_dir", type=str, required=True, help="输入 CIF 文件夹路径")
    parser.add_argument("--output_dir", type=str, required=True, help="输出 PDB 文件夹路径")
    parser.add_argument("--num_processes", type=int, default=None, 
                        help="进程数，默认为CPU数量的80%")

    args = parser.parse_args()
    convert_cif_to_pdb(args.input_dir, args.output_dir, args.num_processes)

