#!/usr/bin/env python3
"""
运行 Hi-C pairs 到 h5ad 转换的便捷脚本
"""

import os
import sys
from pairs_to_h5ad_converter import PairsToH5ADConverter

def main():
    import argparse

    parser = argparse.ArgumentParser(description='Hi-C Pairs to H5AD 转换器 - 便捷运行脚本')
    parser.add_argument('--batch_size', type=int, default=200,
                       help='分批处理大小 (默认: 200, 推荐用于内存优化)')
    parser.add_argument('--use_float64', action='store_true',
                       help='使用 float64 精度 (默认: float32)')
    parser.add_argument('--min_contacts', type=int, default=2,
                       help='最小接触次数阈值 (默认: 2)')
    parser.add_argument('--no_batch', action='store_true',
                       help='禁用分批处理 (需要大量内存)')
    parser.add_argument('--no_chr_file', action='store_true',
                       help='不使用染色体文件，改为自动检测 (可能不完整)')

    args = parser.parse_args()

    # 定义路径
    project_root = "/Volumes/SumSung500/CSU/0_HiRES"
    cell_table_path = os.path.join(project_root, "11pair2h5ad", "cell_table.tsv")
    output_dir = os.path.join(project_root, "11pair2h5ad", "output")
    chr_sizes_file = os.path.join(project_root, "11pair2h5ad", "mm10_chrom_sizes_with_chrY.txt")

    # 检查染色体文件是否存在
    if not args.no_chr_file and not os.path.exists(chr_sizes_file):
        print(f"⚠️  染色体文件不存在: {chr_sizes_file}")
        print("将使用自动检测模式，可能不包含所有染色体")
        chr_sizes_file = None
    elif args.no_chr_file:
        chr_sizes_file = None

    # 根据配置生成输出文件名
    precision = "float64" if args.use_float64 else "float32"
    batch_info = "nobatch" if args.no_batch else f"batch{args.batch_size}"
    chr_info = "ref" if chr_sizes_file else "auto"
    output_filename = f"hic_data_100k_{precision}_{batch_info}_min{args.min_contacts}_{chr_info}.h5ad"
    output_path = os.path.join(output_dir, output_filename)

    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)

    # 检查输入文件
    if not os.path.exists(cell_table_path):
        print(f"错误: 细胞表文件不存在: {cell_table_path}")
        sys.exit(1)

    print("=" * 80)
    print("Hi-C Pairs to H5AD 转换器 - 雌雄混合支持版")
    print("=" * 80)
    print(f"输入文件: {cell_table_path}")
    print(f"输出文件: {output_path}")
    print(f"分辨率: 100K")
    print(f"分批大小: {args.batch_size if not args.no_batch else '不分批'}")
    print(f"数据精度: {precision}")
    print(f"最小接触数: {args.min_contacts}")
    print(f"染色体信息: {'参考文件' if chr_sizes_file else '自动检测'}")
    if chr_sizes_file:
        print(f"染色体文件: {chr_sizes_file}")
    print("=" * 80)

    # 内存使用预估
    if args.no_batch:
        print("⚠️  警告: 不分批模式需要约 32-64GB 内存")
    else:
        estimated_memory = args.batch_size * 150  # MB per batch
        print(f"📊 预估内存使用: ~{estimated_memory}MB (批处理模式)")

    if chr_sizes_file:
        print("✅ 将使用完整的21条染色体 (支持雌雄混合细胞)")
    else:
        print("⚠️  自动检测模式可能遗漏Y染色体数据")

    print("=" * 80)

    # 创建转换器并执行转换
    converter = PairsToH5ADConverter(
        cell_table_path=cell_table_path,
        output_path=output_path,
        resolution=100000,
        batch_size=None if args.no_batch else args.batch_size,
        use_float32=not args.use_float64,
        min_contacts=args.min_contacts,
        chr_sizes_file=chr_sizes_file
    )

    try:
        adata = converter.convert_to_h5ad()
        print("\n🎉 转换成功完成!")
        print(f"📁 数据已保存到: {output_path}")
        print(f"📊 数据维度: {adata.shape}")
        print(f"💾 数据类型: {adata.X.dtype}")

        # 计算文件大小
        file_size = os.path.getsize(output_path) / (1024**3)  # GB
        print(f"📦 文件大小: {file_size:.2f} GB")

    except Exception as e:
        print(f"\n❌ 转换过程中出现错误: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()