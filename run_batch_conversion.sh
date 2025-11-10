#!/bin/sh
#SBATCH --job-name=cif2pdb_batch
#SBATCH --output=logs/cif2pdb_%A_%a.out
#SBATCH --error=logs/cif2pdb_%A_%a.err
#SBATCH --array=25          # 40批次，根据实际文件数调整
#SBATCH --cpus-per-task=1

# 创建日志目录
mkdir -p logs output_pdbs

# 配置参数
INPUT_DIR="/home/megagatlingpea/xiaoyifei/"
OUTPUT_DIR="/home/megagatlingpea/yifei/output_pdbs"
BATCH_SIZE=5000
SAVE_MODE="separate"  # 可选: first, all, separate
SANITY_CHECK="true"

# Python脚本
python3 << 'PYTHON_SCRIPT'
import os
import sys
import glob
import json
from cif2pdb_batch import parse_complete_cif_batch

task_id = int(os.environ.get('SLURM_ARRAY_TASK_ID', 0))
batch_size = int(os.environ.get('BATCH_SIZE', 5000))
input_dir = os.environ.get('INPUT_DIR', '/home/megagatlingpea/xiaoyifei/')
output_dir = os.environ.get('OUTPUT_DIR', './output_pdbs')
save_mode = os.environ.get('SAVE_MODE', 'separate')
sanity_check = os.environ.get('SANITY_CHECK', 'true').lower() == 'true'

# 获取所有CIF文件
all_cifs = sorted(glob.glob(f"{input_dir}/**/*.cif.gz", recursive=True))
print(f"Total CIF files found: {len(all_cifs)}")

# 计算本批次处理的文件范围
start_idx = task_id * batch_size
end_idx = min(start_idx + batch_size, len(all_cifs))
batch_files = all_cifs[start_idx:end_idx]

print(f"Task {task_id}: Processing files {start_idx} to {end_idx-1} ({len(batch_files)} files)")

if not batch_files:
    print(f"No files to process for task {task_id}")
    sys.exit(0)

# 创建批次专属输出目录
batch_output = os.path.join(output_dir, f"batch_{task_id:03d}")

# 执行转换
results = parse_complete_cif_batch(
    batch_files,
    dump_root=batch_output,
    save_mode=save_mode,
    sanity_check=sanity_check
)

# 保存结果日志
results_file = f"logs/results_batch_{task_id:03d}.json"
with open(results_file, 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nBatch {task_id} completed:")
print(f"  Processed: {results['processed']}")
print(f"  Success: {len(results['success'])}")
print(f"  Failed: {len(results['failed'])}")
print(f"Results saved to: {results_file}")
PYTHON_SCRIPT

echo "Task $SLURM_ARRAY_TASK_ID completed at $(date)"
