#!/usr/bin/env python3
"""
测试修复后的 bin 格式插补函数
"""

import os
import sys
import logging

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def test_fixed_imputation():
    """测试修复后的插补函数"""
    
    # 添加路径
    sys.path.insert(0, '/Volumes/SumSung500/CSU/0_HiRES/lwh_schicluster_impute_code/imputation_mac/core')
    sys.path.insert(0, '/Volumes/SumSung500/CSU/0_HiRES/schicluster')
    
    try:
        from impute_bin_format import impute_bin_matrix
        logging.info("✅ 成功导入 impute_bin_format")
    except Exception as e:
        logging.error(f"❌ 导入失败: {e}")
        return False
    
    # 测试文件
    test_file = "/Volumes/SumSung500/CSU/0_HiRES/output/pairs2matrix_output/E70/GasfE703153/GasfE703153_chr1.txt"
    chrom_size_file = "/Volumes/SumSung500/CSU/0_HiRES/mm10_chrom_sizes.txt"
    output_file = "/tmp/test_impute_chr1.npz"
    
    logging.info(f"🧪 测试文件: {test_file}")
    
    if not os.path.exists(test_file):
        logging.error(f"❌ 测试文件不存在: {test_file}")
        return False
    
    if not os.path.exists(chrom_size_file):
        logging.error(f"❌ 染色体大小文件不存在: {chrom_size_file}")
        return False
    
    try:
        success = impute_bin_matrix(
            bin_contact_file=test_file,
            output_path=output_file,
            chrom_size_file=chrom_size_file,
            chrom='chr1',
            resolution=20000,
            use_mps=True
        )
        
        if success:
            logging.info("🎉 测试成功!")
            if os.path.exists(output_file):
                size = os.path.getsize(output_file)
                logging.info(f"📄 输出文件大小: {size:,} bytes")
            return True
        else:
            logging.error("❌ 插补函数返回失败")
            return False
            
    except Exception as e:
        logging.error(f"❌ 测试过程出错: {e}")
        import traceback
        logging.error(traceback.format_exc())
        return False

if __name__ == "__main__":
    logging.info("🚀 开始测试修复后的插补函数")
    test_fixed_imputation()
