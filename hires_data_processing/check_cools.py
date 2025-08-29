import cooler
import numpy as np

cool_path = '/Volumes/SumSung500/CSU/0_HiRES/hires_data_processing/outputs/E70/impute/100K/chunk0/GSM6999841_GasdE701331.cool'

try:
    c = cooler.Cooler(cool_path)
    
    print('=== .cool文件详细分析 ===')
    print(f'文件路径: {cool_path}')
    print(f'分辨率: {c.binsize:,} bp')
    print(f'总bin数: {c.bins().shape[0]:,}')
    print(f'非零像素数: {c.pixels().shape[0]:,}')
    print(f'染色体: {list(c.chromnames)}')
    
    # 检查几个主要染色体的接触情况
    test_chroms = ['chr1', 'chr2', 'chr3']
    total_nonzero = 0
    total_elements = 0
    
    for chrom in test_chroms:
        if chrom in c.chromnames:
            extent = c.extent(chrom)
            chrom_bins = extent[1] - extent[0]
            
            # 取该染色体的前50x50子矩阵进行检查
            test_size = min(50, chrom_bins)
            matrix = c.matrix(balance=False)[extent[0]:extent[0]+test_size, extent[0]:extent[0]+test_size]
            
            nonzero = (matrix > 0).sum()
            total_nonzero += nonzero
            total_elements += matrix.size
            
            print(f'{chrom}: {test_size}x{test_size}矩阵, 非零元素: {nonzero}/{matrix.size} ({100*nonzero/matrix.size:.1f}%)')
            
            # 显示一些非零值
            if nonzero > 0:
                nonzero_vals = matrix[matrix > 0]
                print(f'  非零值范围: {nonzero_vals.min():.4f} - {nonzero_vals.max():.4f}')
                print(f'  前5个非零值: {nonzero_vals[:5]}')
    
    print(f'\\n总体统计: {total_nonzero}/{total_elements} 非零元素 ({100*total_nonzero/total_elements:.1f}%)')
    
    # 结论
    if total_nonzero == 0:
        print('\\n⚠️  结论: 上三角接触对全部为0，可能是空矩阵或插补失败')
    else:
        print(f'\\n✓ 结论: 包含{total_nonzero}个非零接触对，数据正常')
        
except Exception as e:
    print(f'错误: {e}')