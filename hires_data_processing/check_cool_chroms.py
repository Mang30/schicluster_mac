import cooler
import os
import sys
import glob
import pandas as pd

def check_cool_files(directory):
    """
    Checks all .cool files in a directory and reports the number of chromosomes.
    """
    cool_files = sorted(glob.glob(os.path.join(directory, '*.cool')))
    if not cool_files:
        print(f"在目录 {directory} 中没有找到 .cool 文件")
        return

    print(f"找到 {len(cool_files)} 个 .cool 文件。正在检查染色体数量...")
    print("-" * 60)

    results = []
    for cool_file_path in cool_files:
        filename = os.path.basename(cool_file_path)
        try:
            c = cooler.Cooler(cool_file_path)
            num_chroms = len(c.chroms()[:])
            results.append({'File': filename, 'Chromosome Count': num_chroms})
        except Exception as e:
            results.append({'File': filename, 'Chromosome Count': f"错误: {e}"})

    # 使用 pandas 以更美观的表格形式打印结果
    if results:
        df = pd.DataFrame(results)
        print(df.to_string())
    
    print("-" * 60)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("用法: python check_cool_chroms.py <directory_path>")
        sys.exit(1)
    
    target_directory = sys.argv[1]
    if not os.path.isdir(target_directory):
        print(f"错误: 目录不存在 {target_directory}")
        sys.exit(1)
        
    check_cool_files(target_directory)
