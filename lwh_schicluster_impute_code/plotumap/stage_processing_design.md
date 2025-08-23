# Stage分组处理方案设计

## 1. 文件匹配规则

根据分析，scHi-C数据文件和metadata中的Cellname匹配规则如下：

1. scHi-C文件命名格式：`GSMXXXXXX_Cellname.pairs.gz` (例如: `GSM6998595_GasaE751001.pairs.gz`)
2. metadata中的Cellname字段：直接对应文件名中的Cellname部分 (例如: `GasaE751001`)
3. 匹配方法：从文件名中提取`_`之后、`.pairs.gz`之前的字符串作为Cellname进行匹配

## 2. Stage分组策略

根据metadata分析，Stage包含以下分类：
- E70 (557个细胞)
- E75 (1870个细胞)
- E80 (559个细胞)
- E85 (1125个细胞)
- E95 (1108个细胞)
- EX05 (1128个细胞)
- EX15 (1122个细胞)

## 3. 处理流程设计

### 3.1 数据预处理阶段
1. 读取metadata文件，获取所有细胞的Cellname、Stage和Celltype信息
2. 遍历GSE223917_RAW目录下的所有.pairs.gz文件
3. 根据文件名提取Cellname，与metadata进行匹配
4. 按Stage分组存储文件路径和对应的细胞信息

### 3.2 Stage处理阶段
对每个Stage:
1. 读取该Stage下所有细胞的.pairs.gz文件
2. 解析Hi-C接触数据，构建接触矩阵
3. 将接触矩阵转换为h5ad格式
4. 添加细胞类型注释信息
5. 保存Stage特定的h5ad文件

### 3.3 UMAP可视化阶段
1. 对每个Stage的h5ad文件进行UMAP降维分析
2. 使用统一的颜色映射方案显示细胞类型
3. 生成各Stage的UMAP图

## 4. 技术实现细节

### 4.1 h5ad文件格式要求
- 保持与K562_T1_2k_sim.h5ad一致的稀疏矩阵格式
- obs中包含Cellname、Stage、Celltype等元数据
- var中包含基因组bin信息
- layers中包含counts数据

### 4.2 文件处理优化
- 由于数据量大(15793个文件)，需要分批处理
- 可以先处理数据量较小的Stage进行验证(E70或E80)
- 实现断点续传机制，避免重复处理已处理的文件

## 5. 输出文件结构
```
output/
├── stage_E70.h5ad
├── stage_E75.h5ad
├── stage_E80.h5ad
├── stage_E85.h5ad
├── stage_E95.h5ad
├── stage_EX05.h5ad
├── stage_EX15.h5ad
├── umap_E70.png
├── umap_E75.png
├── umap_E80.png
├── umap_E85.png
├── umap_E95.png
├── umap_EX05.png
├── umap_EX15.png
└── color_mapping.json
```