#!/usr/bin/env python3
"""
Hi-C Pairs to H5AD Converter

将细胞 Hi-C pairs 数据转换为 h5ad 格式，使用 100K 分辨率
作者: Claude Code
"""

import os
import sys
import gzip
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix
from collections import defaultdict
import argparse
from tqdm import tqdm
import logging

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class PairsToH5ADConverter:
    def __init__(self, cell_table_path, output_path, resolution=100000,
                 batch_size=None, use_float32=True, min_contacts=1,
                 chr_sizes_file=None):
        """
        初始化转换器

        Args:
            cell_table_path: cell_table.tsv 文件路径
            output_path: 输出 h5ad 文件路径
            resolution: 分辨率，默认 100K
            batch_size: 分批处理大小，None表示不分批
            use_float32: 是否使用float32降低内存使用
            min_contacts: 最小接触次数阈值，低于此值的接触将被过滤
            chr_sizes_file: 染色体大小文件路径，如果提供则使用此文件而非从pairs文件推断
        """
        self.cell_table_path = cell_table_path
        self.output_path = output_path
        self.resolution = resolution
        self.batch_size = batch_size
        self.use_float32 = use_float32
        self.min_contacts = min_contacts
        self.chr_sizes_file = chr_sizes_file
        self.chromosomes = []
        self.chr_lengths = {}
        self.bin_mapping = {}
        self.total_bins = 0

    def load_cell_table(self):
        """加载细胞表"""
        logger.info(f"加载细胞表: {self.cell_table_path}")
        self.cell_df = pd.read_csv(self.cell_table_path, sep='\t', header=None,
                                  names=['cell_id', 'file_path'])
        logger.info(f"发现 {len(self.cell_df)} 个细胞")
        return self.cell_df

    def get_chromosome_info(self, pairs_file):
        """从 pairs 文件头部获取染色体信息"""
        chr_info = {}

        with gzip.open(pairs_file, 'rt') as f:
            for line in f:
                if line.startswith('#chromosome:'):
                    parts = line.strip().split()
                    chr_name = parts[1]
                    chr_length = int(parts[2])
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
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        chr_name = parts[0]
                        chr_length = int(parts[1])
                        chr_info[chr_name] = chr_length

        logger.info(f"从文件加载了 {len(chr_info)} 条染色体信息")
        return chr_info

    def scan_multiple_files_for_chromosomes(self, cell_df, max_files=10):
        """扫描多个文件以获取完整的染色体信息"""
        logger.info("扫描多个文件以获取完整染色体信息...")

        all_chr_info = {}
        files_scanned = 0

        for idx, row in cell_df.iterrows():
            if files_scanned >= max_files:
                break

            pairs_file = row['file_path']
            if not os.path.exists(pairs_file):
                continue

            try:
                chr_info = self.get_chromosome_info(pairs_file)
                # 合并染色体信息，保留最大长度
                for chr_name, chr_length in chr_info.items():
                    if chr_name not in all_chr_info:
                        all_chr_info[chr_name] = chr_length
                    else:
                        all_chr_info[chr_name] = max(all_chr_info[chr_name], chr_length)

                files_scanned += 1
                logger.debug(f"扫描文件 {files_scanned}: {os.path.basename(pairs_file)}")

            except Exception as e:
                logger.warning(f"扫描文件 {pairs_file} 时出错: {e}")
                continue

        logger.info(f"扫描了 {files_scanned} 个文件，发现 {len(all_chr_info)} 条染色体")
        return all_chr_info

    def setup_bins(self, cell_df=None):
        """建立基因组分箱系统"""
        logger.info("设置基因组分箱...")

        # 优先使用提供的染色体大小文件
        if self.chr_sizes_file and os.path.exists(self.chr_sizes_file):
            chr_info = self.load_chromosome_info_from_file(self.chr_sizes_file)
        elif cell_df is not None:
            # 扫描多个文件获取完整染色体信息
            chr_info = self.scan_multiple_files_for_chromosomes(cell_df)
        else:
            # 回退到原有方法（使用第一个文件）
            sample_file = self.cell_df.iloc[0]['file_path']
            if not os.path.exists(sample_file):
                raise FileNotFoundError(f"示例文件不存在: {sample_file}")
            chr_info = self.get_chromosome_info(sample_file)
            logger.warning("使用第一个文件推断染色体信息，可能不完整")

        if not chr_info:
            raise ValueError("无法获取染色体信息")

        # 按染色体名称排序 - 改进的排序逻辑
        def sort_chromosome(chr_name):
            # 移除 'chr' 前缀进行排序
            if chr_name.startswith('chr'):
                suffix = chr_name[3:]
                if suffix.isdigit():
                    return (0, int(suffix))  # 常染色体按数字排序
                elif suffix == 'X':
                    return (1, 0)  # X染色体
                elif suffix == 'Y':
                    return (1, 1)  # Y染色体
                else:
                    return (2, suffix)  # 其他染色体按字母排序
            else:
                return (3, chr_name)  # 无chr前缀的染色体

        sorted_chrs = sorted(chr_info.keys(), key=sort_chromosome)

        self.chromosomes = sorted_chrs
        self.chr_lengths = chr_info

        bin_start = 0
        for chr_name in self.chromosomes:
            chr_length = chr_info[chr_name]
            n_bins = (chr_length + self.resolution - 1) // self.resolution  # 向上取整

            self.bin_mapping[chr_name] = {
                'start_bin': bin_start,
                'n_bins': n_bins,
                'length': chr_length
            }
            bin_start += n_bins

        self.total_bins = bin_start

        # 详细输出染色体信息
        logger.info(f"总共 {self.total_bins} 个分箱，覆盖 {len(self.chromosomes)} 条染色体")
        logger.info(f"染色体列表: {', '.join(self.chromosomes)}")

        # 特别标注性染色体
        sex_chrs = [chr_name for chr_name in self.chromosomes if chr_name in ['chrX', 'chrY']]
        if sex_chrs:
            logger.info(f"包含性染色体: {', '.join(sex_chrs)}")

        if 'chrY' in self.chromosomes:
            logger.info("✓ 支持雄性细胞 (包含chrY)")
        else:
            logger.warning("⚠ 仅支持雌性细胞 (缺少chrY)")

    def pos_to_bin(self, chr_name, pos):
        """将染色体位置转换为分箱索引"""
        if chr_name not in self.bin_mapping:
            return None

        local_bin = pos // self.resolution
        global_bin = self.bin_mapping[chr_name]['start_bin'] + local_bin

        # 确保不超出该染色体的范围
        if local_bin >= self.bin_mapping[chr_name]['n_bins']:
            return None

        return global_bin

    def process_pairs_file(self, pairs_file):
        """处理单个 pairs 文件，返回接触矩阵"""
        logger.debug(f"处理文件: {pairs_file}")

        contacts = defaultdict(int)

        try:
            with gzip.open(pairs_file, 'rt') as f:
                for line in f:
                    if line.startswith('#'):
                        continue

                    parts = line.strip().split('\t')
                    if len(parts) < 4:
                        continue

                    chr1, pos1, chr2, pos2 = parts[1], int(parts[2]), parts[3], int(parts[4])

                    bin1 = self.pos_to_bin(chr1, pos1)
                    bin2 = self.pos_to_bin(chr2, pos2)

                    if bin1 is not None and bin2 is not None:
                        # 确保 bin1 <= bin2 (上三角矩阵)
                        if bin1 > bin2:
                            bin1, bin2 = bin2, bin1
                        contacts[(bin1, bin2)] += 1

        except Exception as e:
            logger.error(f"处理文件 {pairs_file} 时出错: {e}")
            return None

        # 转换为稀疏矩阵
        if not contacts:
            logger.warning(f"文件 {pairs_file} 中没有找到有效的接触数据")
            triu_indices = np.triu_indices(self.total_bins)
            n_features = len(triu_indices[0])
            dtype = np.float32 if self.use_float32 else np.float64
            return np.zeros(n_features, dtype=dtype)

        # 过滤低接触次数
        filtered_contacts = {k: v for k, v in contacts.items() if v >= self.min_contacts}

        if not filtered_contacts:
            logger.warning(f"文件 {pairs_file} 过滤后没有有效的接触数据")
            triu_indices = np.triu_indices(self.total_bins)
            n_features = len(triu_indices[0])
            dtype = np.float32 if self.use_float32 else np.float64
            return np.zeros(n_features, dtype=dtype)

        # 直接构建上三角向量，节省内存
        triu_indices = np.triu_indices(self.total_bins)
        n_features = len(triu_indices[0])
        dtype = np.float32 if self.use_float32 else np.float64
        contact_vector = np.zeros(n_features, dtype=dtype)

        # 创建索引映射字典（从(i,j)到向量索引）
        index_map = {}
        for idx, (i, j) in enumerate(zip(triu_indices[0], triu_indices[1])):
            index_map[(i, j)] = idx

        # 填充接触数据
        for (bin1, bin2), count in filtered_contacts.items():
            if (bin1, bin2) in index_map:
                contact_vector[index_map[(bin1, bin2)]] = count

        return contact_vector

    def process_batch(self, cell_batch):
        """处理一批细胞"""
        batch_matrices = []
        batch_cells = []

        for idx, row in cell_batch.iterrows():
            cell_id = row['cell_id']
            pairs_file = row['file_path']

            if not os.path.exists(pairs_file):
                logger.warning(f"文件不存在，跳过: {pairs_file}")
                continue

            contact_vector = self.process_pairs_file(pairs_file)

            if contact_vector is not None:
                batch_matrices.append(contact_vector)
                batch_cells.append(cell_id)

        return batch_matrices, batch_cells

    def convert_to_h5ad(self):
        """主转换函数，支持分批处理"""
        logger.info("开始转换 pairs 文件到 h5ad 格式...")

        # 加载细胞表
        cell_df = self.load_cell_table()

        # 设置基因组分箱系统 - 传入cell_df以支持多文件扫描
        self.setup_bins(cell_df)

        # 确定是否使用分批处理
        if self.batch_size is None:
            return self._convert_all_at_once(cell_df)
        else:
            return self._convert_in_batches(cell_df)

    def _convert_all_at_once(self, cell_df):
        """一次性处理所有细胞"""
        logger.info("使用一次性处理模式...")

        # 处理所有细胞
        cell_matrices = []
        valid_cells = []

        for idx, row in tqdm(cell_df.iterrows(), total=len(cell_df), desc="处理细胞"):
            cell_id = row['cell_id']
            pairs_file = row['file_path']

            if not os.path.exists(pairs_file):
                logger.warning(f"文件不存在，跳过: {pairs_file}")
                continue

            contact_vector = self.process_pairs_file(pairs_file)

            if contact_vector is not None:
                cell_matrices.append(contact_vector)
                valid_cells.append(cell_id)

            # 每处理100个细胞输出一次进度
            if (idx + 1) % 100 == 0:
                logger.info(f"已处理 {idx + 1}/{len(cell_df)} 个细胞")

        if not cell_matrices:
            raise ValueError("没有成功处理任何细胞数据")

        logger.info(f"成功处理 {len(cell_matrices)} 个细胞")
        return self._create_anndata(cell_matrices, valid_cells)

    def _convert_in_batches(self, cell_df):
        """分批处理细胞"""
        logger.info(f"使用分批处理模式，批大小: {self.batch_size}")

        all_matrices = []
        all_cells = []
        n_batches = (len(cell_df) + self.batch_size - 1) // self.batch_size

        for batch_idx in range(n_batches):
            start_idx = batch_idx * self.batch_size
            end_idx = min((batch_idx + 1) * self.batch_size, len(cell_df))
            batch_df = cell_df.iloc[start_idx:end_idx]

            logger.info(f"处理批次 {batch_idx + 1}/{n_batches} ({len(batch_df)} 个细胞)")

            batch_matrices, batch_cells = self.process_batch(batch_df)

            if batch_matrices:
                all_matrices.extend(batch_matrices)
                all_cells.extend(batch_cells)

            # 显示进度
            logger.info(f"批次 {batch_idx + 1} 完成，本批处理 {len(batch_matrices)} 个细胞")

            # 激进内存清理
            del batch_matrices, batch_cells
            import gc
            gc.collect()

            # 强制内存释放
            import psutil
            process = psutil.Process()
            logger.debug(f"批次 {batch_idx + 1} 完成后内存使用: {process.memory_info().rss / 1024**3:.2f} GB")

        if not all_matrices:
            raise ValueError("没有成功处理任何细胞数据")

        logger.info(f"所有批次处理完成，共处理 {len(all_matrices)} 个细胞")
        return self._create_anndata(all_matrices, all_cells)

    def _create_anndata(self, cell_matrices, valid_cells):
        """创建 AnnData 对象"""
        logger.info("创建 AnnData 对象...")

        # 转换为numpy数组
        dtype = np.float32 if self.use_float32 else np.float64
        X = np.vstack(cell_matrices).astype(dtype)

        # 创建观测值注释（细胞）
        obs = pd.DataFrame(index=valid_cells)
        obs['cell_type'] = 'unknown'

        # 创建变量注释（基因组区间对）
        n_features = len(cell_matrices[0])
        var = pd.DataFrame(index=[f'contact_{i}' for i in range(n_features)])
        var['feature_type'] = 'chromatin_contact'

        # 创建 AnnData 对象
        adata = ad.AnnData(X=X, obs=obs, var=var)

        # 添加元数据
        adata.uns['resolution'] = self.resolution
        adata.uns['chromosomes'] = self.chromosomes
        adata.uns['chr_lengths'] = self.chr_lengths
        adata.uns['bin_mapping'] = self.bin_mapping
        adata.uns['total_bins'] = self.total_bins
        adata.uns['min_contacts'] = self.min_contacts
        adata.uns['use_float32'] = self.use_float32

        # 保存
        logger.info(f"保存 h5ad 文件: {self.output_path}")
        adata.write_h5ad(self.output_path)

        logger.info(f"转换完成! 数据形状: {adata.shape}")
        logger.info(f"细胞数量: {adata.n_obs}")
        logger.info(f"特征数量: {adata.n_vars}")
        logger.info(f"数据类型: {adata.X.dtype}")

        return adata

def main():
    parser = argparse.ArgumentParser(description='将 Hi-C pairs 数据转换为 h5ad 格式')
    parser.add_argument('--cell_table', required=True, help='cell_table.tsv 文件路径')
    parser.add_argument('--output', required=True, help='输出 h5ad 文件路径')
    parser.add_argument('--resolution', type=int, default=100000, help='分辨率 (默认: 100000)')
    parser.add_argument('--batch_size', type=int, default=None, help='分批处理大小 (默认: None, 不分批)')
    parser.add_argument('--use_float64', action='store_true', help='使用 float64 精度 (默认: float32)')
    parser.add_argument('--min_contacts', type=int, default=1, help='最小接触次数阈值 (默认: 1)')
    parser.add_argument('--chr_sizes', type=str, help='染色体大小文件路径 (推荐使用以确保完整性)')
    parser.add_argument('--verbose', action='store_true', help='详细输出')

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # 检查输入文件
    if not os.path.exists(args.cell_table):
        print(f"错误: 细胞表文件不存在: {args.cell_table}")
        sys.exit(1)

    # 检查染色体大小文件
    if args.chr_sizes and not os.path.exists(args.chr_sizes):
        print(f"错误: 染色体大小文件不存在: {args.chr_sizes}")
        sys.exit(1)

    # 创建输出目录
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 输出配置信息
    print("=" * 60)
    print("Hi-C Pairs to H5AD 转换器配置")
    print("=" * 60)
    print(f"分辨率: {args.resolution:,}")
    print(f"分批大小: {args.batch_size if args.batch_size else '不分批'}")
    print(f"数据精度: {'float64' if args.use_float64 else 'float32'}")
    print(f"最小接触数: {args.min_contacts}")
    print(f"染色体信息: {args.chr_sizes if args.chr_sizes else '自动检测'}")
    print("=" * 60)

    # 执行转换
    converter = PairsToH5ADConverter(
        cell_table_path=args.cell_table,
        output_path=args.output,
        resolution=args.resolution,
        batch_size=args.batch_size,
        use_float32=not args.use_float64,
        min_contacts=args.min_contacts,
        chr_sizes_file=args.chr_sizes
    )
    converter.convert_to_h5ad()

if __name__ == '__main__':
    main()