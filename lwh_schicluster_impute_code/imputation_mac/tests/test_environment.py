#!/usr/bin/env python3
"""
M4 MacBook 环境测试脚本
验证MPS加速、数据路径和依赖项
"""

import os
import sys
import logging

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def test_mps_availability():
    """测试MPS加速可用性"""
    try:
        import torch
        logging.info(f"✅ PyTorch 版本: {torch.__version__}")
        
        if torch.backends.mps.is_available():
            logging.info("✅ MPS 加速可用")
            
            # 测试MPS设备创建
            device = torch.device('mps')
            test_tensor = torch.randn(10, 10, device=device)
            result = torch.mm(test_tensor, test_tensor.t())
            logging.info(f"✅ MPS 设备测试成功，输出形状: {result.shape}")
            return True
        else:
            logging.warning("⚠️  MPS 不可用，将使用 CPU")
            return False
            
    except ImportError as e:
        logging.error(f"❌ PyTorch 导入失败: {e}")
        return False
    except Exception as e:
        logging.error(f"❌ MPS 测试失败: {e}")
        return False

def test_schicluster_import():
    """测试schicluster导入"""
    try:
        sys.path.insert(0, '/Volumes/SumSung500/CSU/0_HiRES/schicluster')
        from schicluster.impute.impute_chromosome import impute_chromosome, MPS_AVAILABLE
        logging.info("✅ schicluster 导入成功")
        logging.info(f"✅ schicluster MPS 支持: {MPS_AVAILABLE}")
        return True
    except ImportError as e:
        logging.error(f"❌ schicluster 导入失败: {e}")
        return False

def test_data_paths():
    """测试数据路径"""
    paths_to_check = {
        "输入数据": "/Volumes/SumSung500/CSU/0_HiRES/output/pairs2matrix_output",
        "输出目录": "/Volumes/SumSung500/CSU/0_HiRES/output/imputed_matrices_by_stage",
        "染色体长度": "/Volumes/SumSung500/CSU/0_HiRES/mm10_chrom_sizes.txt",
        "脚本目录": "/Volumes/SumSung500/CSU/0_HiRES/lwh_schicluster_impute_code/imputation_mac"
    }
    
    all_good = True
    for name, path in paths_to_check.items():
        if os.path.exists(path):
            logging.info(f"✅ {name}: {path}")
        else:
            logging.error(f"❌ {name}不存在: {path}")
            all_good = False
    
    return all_good

def test_sample_data():
    """测试样本数据"""
    try:
        input_dir = "/Volumes/SumSung500/CSU/0_HiRES/output/pairs2matrix_output"
        stages = ["E70", "E75", "E80", "E85", "E95", "EX05", "EX15"]
        
        total_cells = 0
        for stage in stages:
            stage_path = os.path.join(input_dir, stage)
            if os.path.exists(stage_path):
                cells = [d for d in os.listdir(stage_path) 
                        if (d.startswith("Gas") or d.startswith("Org")) and os.path.isdir(os.path.join(stage_path, d))]
                cell_count = len(cells)
                total_cells += cell_count
                logging.info(f"✅ {stage}: {cell_count} 个细胞")
                
                # 检查第一个细胞的数据
                if cells:
                    first_cell = cells[0]
                    cell_path = os.path.join(stage_path, first_cell)
                    chr_files = [f for f in os.listdir(cell_path) if f.endswith('.txt')]
                    logging.info(f"   样本细胞 {first_cell}: {len(chr_files)} 个染色体文件")
            else:
                logging.warning(f"⚠️  {stage} 目录不存在")
        
        logging.info(f"✅ 总计: {total_cells} 个细胞")
        return total_cells > 0
        
    except Exception as e:
        logging.error(f"❌ 样本数据检查失败: {e}")
        return False

def test_memory_info():
    """检查内存信息"""
    try:
        import psutil
        memory = psutil.virtual_memory()
        logging.info(f"💻 系统内存: {memory.total / (1024**3):.1f}GB 总计")
        logging.info(f"💻 可用内存: {memory.available / (1024**3):.1f}GB")
        logging.info(f"💻 内存使用: {memory.percent:.1f}%")
        return True
    except ImportError:
        logging.warning("⚠️  psutil 未安装，无法检查内存信息")
        return True
    except Exception as e:
        logging.error(f"❌ 内存信息检查失败: {e}")
        return False

def main():
    """主测试函数"""
    logging.info("🚀 开始 M4 MacBook MPS 插补环境测试")
    logging.info("=" * 60)
    
    tests = [
        ("MPS 加速", test_mps_availability),
        ("schicluster 导入", test_schicluster_import), 
        ("数据路径", test_data_paths),
        ("样本数据", test_sample_data),
        ("内存信息", test_memory_info)
    ]
    
    passed = 0
    failed = 0
    
    for test_name, test_func in tests:
        logging.info(f"\n🔍 测试: {test_name}")
        try:
            if test_func():
                passed += 1
                logging.info(f"✅ {test_name} 通过")
            else:
                failed += 1
                logging.error(f"❌ {test_name} 失败")
        except Exception as e:
            failed += 1
            logging.error(f"❌ {test_name} 异常: {e}")
    
    logging.info("=" * 60)
    logging.info(f"📊 测试结果: {passed} 通过, {failed} 失败")
    
    if failed == 0:
        logging.info("🎉 所有测试通过! 环境已就绪")
        logging.info("📋 下一步: 运行测试插补")
        logging.info("   python batch_mps_imputation_v2.py --stages E70 --max-cells 1")
        return 0
    else:
        logging.error("❌ 部分测试失败，请修复后重试")
        return 1

if __name__ == "__main__":
    exit(main())