这个文件夹下的项目是用于将处理后的各个stage 下的细胞contact decay profile 数据【json 格式】进行h5ad 格式的转换。

【JSON文件结构】
{
  "distances": [2, 3, 4, ...],              // 整数距离bins
  "decay_values": [0.0314, 0.0240, ...],    // Contact decay值
  "distance_kb": [200.0, 300.0, ...],       // 距离的kb表示
  "total_contacts": 10,                      // 总contact数
  "resolution": 100000,                      // 分辨率
  "power_law_slope": -1.045                  // 幂律斜率
}
【文件处理规则】
1. **过滤文件**: 忽略 `analysis_summary.json` 文件，只处理 `*_decay_profile.json` 文件
2. **细胞ID提取**: 从文件名中提取细胞ID（去除 `_decay_profile.json` 后缀）
3. **数据对齐**: 确保所有细胞的距离bins一致，处理长度不同的情况
4. **缺失值处理**: 对于不同长度的decay_values，使用NaN填充或截断到最小长度

【输入数据】要进行格式转换的 stage数据的路径是：
1 `3_create_contact_decay_profile/outputs/stage_E70`
2 `3_create_contact_decay_profile/outputs/stage_E75`
3 `3_create_contact_decay_profile/outputs/stage_E80`
4 `3_create_contact_decay_profile/outputs/stage_E85`
5 `3_create_contact_decay_profile/outputs/stage_E95`
6 `3_create_contact_decay_profile/outputs/stage_EX05`
7 `3_create_contact_decay_profile/outputs/stage_EX15`

【输出数据】保存在 '4_contact_decay_profile_2_h5ad/json2h5ad/output',你需要将同一个 stage 下各个 schic 的contact decay profile 转换为一个 h5ad 文件，并按照stage 分开存放

【创建分开启动的脚本】你需要创建用于将 schic 的contact decay profile 转换为 h5ad 格式的核心转换代码。
同时，你需要创建 7 个启动脚本，分别用于 7 个 stage 下的数据转换。
【环境说明】
使用此环境 micromamba activate schicluster
MAC 下 使用 claude code 运行 micromamba 下环境中的 python 使用此路径： /Users/wuhaoliu/mamba/envs/schicluster/bin/python