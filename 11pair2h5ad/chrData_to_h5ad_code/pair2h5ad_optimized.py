#!/usr/bin/env python3
"""
Hi-C Pairs to H5AD Converter (Optimized Version)

将细胞 Hi-C pairs 数据转换为 h5ad 格式，使用稀疏矩阵和并行处理进行加速。
作者: Claude Code (根据用户代码优化)
"""

import os
import sys
import gzip
import pandas as pd
import numpy as np
import anndata as ad
from scipy.sparse import csr_matrix, vstack
from collections import defaultdict
import argparse
from tqdm import tqdm
import logging
from multiprocessing import Pool
import psutil # 用于获取CPU核心数

# --- 全局函数，用于并行处理 ---
# OPTIMIZATION: 将核心处理逻辑移出类，使其成为一个独立的顶层函数，
# 这样 multiprocessing 就可以轻松地“pickle”(序列化)它。

def _process_single_file_sparse(args):
    """
    处理单个 pairs 文件，返回稀疏接触向量和细胞ID。
    这是一个独立的函数，专为并行化设计。
    """
    # 1. 解包参数
    file_path, resolution, bin_mapping, coord_to_feature_idx, min_contacts, use_float32, n_features = args
    
    cell_id = os.path.basename(file_path).split('.')[0] # 示例：从文件名获取cell_id
    contacts = defaultdict(int)

    try:
        with gzip.open(file_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue
                
                chr1, pos1, chr2, pos2 = parts[1], int(parts[2]), parts[3], int(parts[4])

                # 辅助函数：将位置转换为全局 bin 索引
                def pos_to_bin(chr_name, pos):
                    if chr_name not in bin_mapping: return None
                    local_bin = pos // resolution
                    if local_bin >= bin_mapping[chr_name]['n_bins']: return None
                    return bin_mapping[chr_name]['start_bin'] + local_bin

                bin1 = pos_to_bin(chr1, pos1)
                bin2 = pos_to_bin(chr2, pos2)

                if bin1 is not None and bin2 is not None:
                    if bin1 > bin2: bin1, bin2 = bin2, bin1
                    contacts[(bin1, bin2)] += 1
    except Exception as e:
        # 在多进程中，打印错误而不是使用logger
        print(f"警告: 处理文件 {file_path} 时出错: {e}", file=sys.stderr)
        return None

    # 2. 过滤低接触并准备稀疏矩阵数据
    if not contacts: return None
    
    data, indices = [], []
    for (b1, b2), count in contacts.items():
        if count >= min_contacts:
            feature_idx = coord_to_feature_idx.get((b1, b2))
            if feature_idx is not None:
                data.append(count)
                indices.append(feature_idx)

    if not data: return None

    # 3. 创建一个 1 x n_features 的稀疏矩阵 (CSR格式)
    dtype = np.float32 if use_float32 else np.float64
    contact_matrix_sparse = csr_matrix((data, indices, [0, len(data)]),
                                       shape=(1, n_features),
                                       dtype=dtype)
    
    return contact_matrix_sparse, cell_id

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class PairsToH5ADConverter:
    def __init__(self, cell_table_path, output_path, resolution=100000,
                 use_float32=True, min_contacts=1, chr_sizes_file=None, workers=1):
        """
        初始化转换器
        """
        self.cell_table_path = cell_table_path
        self.output_path = output_path
        self.resolution = resolution
        self.use_float32 = use_float32
        self.min_contacts = min_contacts
        self.chr_sizes_file = chr_sizes_file
        self.workers = workers # OPTIMIZATION: 新增并行工作进程数参数
        
        self.chromosomes = []
        self.chr_lengths = {}
        self.bin_mapping = {}
        self.total_bins = 0
        
        # OPTIMIZATION: 用于存储预计算结果的属性
        self.coord_to_feature_idx = {}
        self.triu_indices = None
        self.n_features = 0

    def load_cell_table(self):
        """加载细胞表"""
        logger.info(f"加载细胞表: {self.cell_table_path}")
        self.cell_df = pd.read_csv(self.cell_table_path, sep='\t', header=None,
                                  names=['cell_id', 'file_path'])
        logger.info(f"发现 {len(self.cell_df)} 个细胞")
        return self.cell_df

    # ... get_chromosome_info, load_chromosome_info_from_file, scan_multiple_files_for_chromosomes 方法保持不变 ...
    def get_chromosome_info(self, pairs_file):
        """从 pairs 文件头部获取染色体信息"""
        chr_info = {}
        with gzip.open(pairs_file, 'rt') as f:
            for line in f:
                if line.startswith('#chromosome:'):
                    parts = line.strip().split()
                    chr_name, chr_length = parts[1], int(parts[2])
                    chr_info[chr_name] = chr_length
                elif not line.startswith('#'):
                    break
        return chr_info

    def load_chromosome_info_from_file(self, chr_sizes_file):
        """从染色体大小文件加载染色体信息"""
        logger.info(f"从文件加载染色体信息: {chr_sizes_file}")
        chr_info = {}
        with open(chr_sizes_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    chr_info[parts[0]] = int(parts[1])
        return chr_info

    def scan_multiple_files_for_chromosomes(self, cell_df, max_files=10):
        """扫描多个文件以获取完整的染色体信息"""
        all_chr_info = {}
        for _, row in cell_df.head(max_files).iterrows():
            if os.path.exists(row['file_path']):
                try:
                    chr_info = self.get_chromosome_info(row['file_path'])
                    for chr_name, chr_length in chr_info.items():
                        all_chr_info[chr_name] = max(all_chr_info.get(chr_name, 0), chr_length)
                except Exception as e:
                    logger.warning(f"扫描文件 {row['file_path']} 时出错: {e}")
        return all_chr_info

    def setup_bins(self, cell_df=None):
        """建立基因组分箱系统并预计算映射"""
        logger.info("设置基因组分箱...")
        if self.chr_sizes_file and os.path.exists(self.chr_sizes_file):
            chr_info = self.load_chromosome_info_from_file(self.chr_sizes_file)
        elif cell_df is not None:
            chr_info = self.scan_multiple_files_for_chromosomes(cell_df)
        else:
            raise ValueError("无法确定染色体信息")
        
        # 染色体排序逻辑
        def sort_key(chr_name):
            suffix = chr_name.lstrip('chr')
            if suffix.isdigit(): return (0, int(suffix))
            if suffix == 'X': return (1, 0)
            if suffix == 'Y': return (1, 1)
            return (2, suffix)
        
        self.chromosomes = sorted(chr_info.keys(), key=sort_key)
        self.chr_lengths = chr_info

        bin_start = 0
        for chr_name in self.chromosomes:
            chr_length = chr_info[chr_name]
            n_bins = (chr_length + self.resolution - 1) // self.resolution
            self.bin_mapping[chr_name] = {'start_bin': bin_start, 'n_bins': n_bins}
            bin_start += n_bins
        self.total_bins = bin_start

        logger.info(f"总共 {self.total_bins} 个分箱，覆盖 {len(self.chromosomes)} 条染色体")

        # --- OPTIMIZATION: 预计算坐标到一维特征索引的映射 ---
        logger.info("预计算坐标到特征索引的映射...")
        rows, cols = np.triu_indices(self.total_bins)
        self.triu_indices = (rows, cols)
        self.n_features = len(rows)
        self.coord_to_feature_idx = {(i, j): idx for idx, (i, j) in enumerate(zip(rows, cols))}
        logger.info(f"上三角矩阵共有 {self.n_features:,} 个特征")

    def convert_to_h5ad(self):
        """主转换函数"""
        logger.info("开始转换 pairs 文件到 h5ad 格式...")
        cell_df = self.load_cell_table()
        self.setup_bins(cell_df)

        # 准备传递给每个进程的参数
        tasks = []
        for _, row in cell_df.iterrows():
            if os.path.exists(row['file_path']):
                tasks.append((
                    row['file_path'], self.resolution, self.bin_mapping,
                    self.coord_to_feature_idx, self.min_contacts,
                    self.use_float32, self.n_features
                ))

        # --- OPTIMIZATION: 根据 workers 数量选择并行或串行处理 ---
        if self.workers > 1:
            logger.info(f"使用 {self.workers} 个工作进程进行并行处理...")
            results = []
            with Pool(self.workers) as pool:
                # imap_unordered 效率更高，结果返回后立即处理
                for res in tqdm(pool.imap_unordered(_process_single_file_sparse, tasks), total=len(tasks), desc="并行处理细胞"):
                    if res is not None:
                        results.append(res)
        else:
            logger.info("使用单进程串行处理...")
            results = [res for res in tqdm(map(_process_single_file_sparse, tasks), total=len(tasks), desc="串行处理细胞") if res is not None]

        if not results:
            raise ValueError("没有成功处理任何细胞数据")
        
        logger.info(f"成功处理 {len(results)} / {len(cell_df)} 个细胞")
        
        # 按原始 cell_id 顺序对结果进行排序，确保 AnnData 行序正确
        cell_order_map = {cell_id: i for i, cell_id in enumerate(cell_df['cell_id'])}
        results.sort(key=lambda x: cell_order_map.get(x[1], -1))

        cell_matrices_sparse = [res[0] for res in results]
        valid_cells = [res[1] for res in results]

        self._create_anndata(cell_matrices_sparse, valid_cells)

    def _create_anndata(self, cell_matrices_sparse, valid_cells):
        """从稀疏矩阵列表创建并保存 AnnData 对象"""
        logger.info("创建 AnnData 对象...")

        # --- OPTIMIZATION: 使用 scipy.sparse.vstack 合并稀疏矩阵 ---
        X_sparse = vstack(cell_matrices_sparse, format='csr')

        obs = pd.DataFrame(index=valid_cells)
        obs['cell_type'] = 'unknown' # 可根据 cell_table 文件内容填充

        # --- 增强：创建包含详细坐标信息的 var DataFrame ---
        logger.info("创建特征(var)注释...")
        var = pd.DataFrame(index=[f'contact_{i}' for i in range(self.n_features)])
        var['bin1_id'] = self.triu_indices[0]
        var['bin2_id'] = self.triu_indices[1]
        
        adata = ad.AnnData(X=X_sparse, obs=obs, var=var, dtype=X_sparse.dtype)

        adata.uns['resolution'] = self.resolution
        adata.uns['chromosomes'] = self.chromosomes
        adata.uns['chr_lengths'] = self.chr_lengths
        adata.uns['bin_mapping'] = self.bin_mapping
        adata.uns['total_bins'] = self.total_bins
        
        logger.info(f"保存 h5ad 文件: {self.output_path}")
        adata.write_h5ad(self.output_path)

        logger.info("=" * 30 + " 转换完成! " + "=" * 30)
        logger.info(f"数据形状: {adata.shape}")
        logger.info(f"数据类型: {adata.X.dtype}")
        logger.info(f"稀疏度: {adata.X.nnz / (adata.shape[0] * adata.shape[1]) * 100:.4f}%")


def main():
    parser = argparse.ArgumentParser(description='将 Hi-C pairs 数据高效转换为 h5ad 格式')
    parser.add_argument('--cell_table', required=True, help='cell_table.tsv 文件路径')
    parser.add_argument('--output', required=True, help='输出 h5ad 文件路径')
    parser.add_argument('--resolution', type=int, default=100000, help='分辨率 (默认: 100000)')
    parser.add_argument('--chr_sizes', type=str, help='染色体大小文件路径 (强烈推荐)')
    parser.add_argument('--min_contacts', type=int, default=1, help='最小接触次数阈值 (默认: 1)')
    # OPTIMIZATION: 新增 workers 参数
    parser.add_argument('--workers', type=int, default=1, help='用于处理文件的并行工作进程数。默认为1(串行)。设为0则使用所有可用CPU核心。')
    parser.add_argument('--use_float64', action='store_true', help='使用 float64 精度 (默认: float32)')
    
    args = parser.parse_args()

    # --- 增强的参数处理 ---
    if not os.path.exists(args.cell_table):
        logger.error(f"错误: 细胞表文件不存在: {args.cell_table}")
        sys.exit(1)
    if args.chr_sizes and not os.path.exists(args.chr_sizes):
        logger.error(f"错误: 染色体大小文件不存在: {args.chr_sizes}")
        sys.exit(1)
        
    num_workers = args.workers
    if num_workers == 0:
        num_workers = psutil.cpu_count(logical=False) # 使用物理核心数以获得最佳性能
    
    output_dir = os.path.dirname(args.output)
    if output_dir: os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 60)
    print("Hi-C Pairs to H5AD 高性能转换器")
    print("=" * 60)
    print(f"并行工作进程数: {num_workers}")
    # ... 其他配置信息 ...
    print("=" * 60)
    
    converter = PairsToH5ADConverter(
        cell_table_path=args.cell_table,
        output_path=args.output,
        resolution=args.resolution,
        use_float32=not args.use_float64,
        min_contacts=args.min_contacts,
        chr_sizes_file=args.chr_sizes,
        workers=num_workers
    )
    converter.convert_to_h5ad()

if __name__ == '__main__':
    main()