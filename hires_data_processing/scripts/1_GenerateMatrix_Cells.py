import os
import subprocess

# 设置输入和输出目录，染色体大小文件，分辨率等
contact_file_dir = './contact_files'
raw_dir = './raw_matrices'
chromosome_size_file = './chromosome_size_file.txt'
resolution = 1000000  # 1Mb分辨率

# 遍历所有接触文件
for contact_file in os.listdir(contact_file_dir):
    if contact_file.endswith('_contact.txt'):  # 确保是以 "_contact.txt" 结尾的文件
        cell_id = contact_file.split('_contact.txt')[0]

        # 构建命令
        command = [
            'hicluster', 'generatematrix-cell',
            '--infile', os.path.join(contact_file_dir, contact_file),
            '--outdir', raw_dir,
            '--chrom_file', chromosome_size_file,
            '--res', str(resolution),
            '--cell', cell_id,
            '--chr1', '1', '--pos1', '2', '--chr2', '3', '--pos2', '4'  # 假设是这些列
        ]

        # 执行命令
        subprocess.run(command)
