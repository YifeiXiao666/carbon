#!/bin/bash

#SBATCH -p com 指定队列名称
#SBATCH -J cif2pdb 指定作业名称
#SBATCH -N 2 指定要提交的节点数量
#SBATCH -n 1 指定要提交的总核数
#SBATCH -o test.o 指定标准输出文件名
#SBATCH -e test.e 指定错误输出文件名

source activate RAPiGen

echo "Job start！"
python3 /home/megagatlingpea/yifei/cif2pdb_batch.py
echo "Job end！"
