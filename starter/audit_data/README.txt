真实数据来源审计输入（不是合成示例）

GSE245552_sample_manifest.csv 可直接输入 validate_manifest.py。
39 GSM / 16 独立患者 / 13 PT–LM 配对 / 1 四组织患者；16 人全部 development。
同一患者的多个病灶继续共享 donor_uid。lesion_uid 按原始 specimen 编号命名，
包括正常邻近组织；不声称每个正常组织是恶性病灶。

组织由 GEO source_name 核实，不能从 T1/T2 的数字自动猜测：
Pt08_T1/T2 均为 PT，Pt08_T3/T4 为 LM；Pt16_T1 是 LM。

原始临床补充 Table S1 DOI 10.1158/0008-5472.25607328.v1 给出 MSS、
部位、年龄、性别和分期；未给治疗，treatment 保留 unknown。
原始 S1 共 18 人，Pt19/Pt20 不在本 GEO 单细胞矩阵；未硬塞进 manifest。

本轮已完整下载并读完39份矩阵，共1,056,396,888压缩字节、347,049,770条MTX记录。
39份均通过非负整数/feature和barcode维度/坐标唯一性检查，counts_values_verified
均为true；GSE245552_all_matrix_audit.csv/json逐文件保存SHA256及输入QC分位数。
大MTX核查后删除，提供audit_GSE245552_all_matrices.py可按原始URL重跑（约1GB）。

关键限制：33个PT/LM只有29569基因，6个NC/NL有29922基因；肿瘤缺353条，
包括全部MT-基因及很多RPL/RPS/MRPL/MRPS/其他线粒体相关基因。
这些是仍保留整数计数的已缩减基因矩阵，不是完整Cell Ranger原输出。
不能在肿瘤上重新计算percent.mt并将0当好质量；不能做缺失线粒体/核糖体程序比较。
程序评分须按冻结的实际可测基因集合；缺失feature不是生物学0。
NB资格只表示这些未取log的计数可进入相应模型，不表示全基因组、细胞QC、
谱系注释或患者级样本量已经充分。剩余细胞质量和作者过滤链仍须按方案审查。

完整临床属性及文件 URL 保留于 GSE245552_sample_manifest_audited.csv，
GSE245552_patient_clinical_S1.csv 保留原临床表 18 人与 scRNA 是否入组标记。

运行（从 CRLM_starter 目录）：
python validate_manifest.py --manifest audit_data/GSE245552_sample_manifest.csv --outdir audit_data/validated

重跑完整矩阵核验（约1GB下载；默认成功后删除大矩阵）：
python audit_data/audit_GSE245552_all_matrices.py --force


默认离线核心核验（只用本包证据；不联网、不下载矩阵）：
python audit_data/reproduce_core_audit.py

此脚本重新读取完整 metadata 身份列子集、六个164表达表头、两份原245feature文件、
完整正常图谱obs选择列、原临床S1，并重聚合39份历史完整扫描记录。
结果明确标为历史完整矩阵检查记录的再聚合，不会声称本次又扫描了1GB矩阵。
metadata身份子集保留原文件每一行，只留barcode/patient_code/tissue；未重新命名
或筛选患者。identity_subset_provenance.json记录原下载文件SHA及子集SHA。
正常图谱obs只含原标签，没有表达，不代表重新注释或生物效应分析。

需要从GEO重新核验全部245矩阵时，明确执行：
python audit_data/audit_GSE245552_all_matrices.py --force
约1GB下载，默认核验后删除大MTX。只调用脚本而不给--force不会偷偷下载缺失文件。
附带的其他网络探针是可选重跑工具，默认离线入口不会调用它们。
其中probe_normal_atlas_metadata.py需要额外安装h5py并使用网络；它不属于默认
离线核心入口。默认reproduce_core_audit.py只依赖Python标准库及包内脚本。

原始GSE164522全metadata没有重复装入包；若重跑在线六表头审计，先把官方
GSE164522_CRLM_metadata.csv.gz存到audit_data/geo。身份列子集不能冒充完整原文件。
可用URL及原始数据条款见 ../evidence/data_sources_D01_D10.txt。
GSE225857补充表已取得，六人具体治疗/部位已核，但额外s0816身份仍未解决。
六个临床代码都不能凭日期格式自动合并，MSS/MMR在S2-S4没有列。

245原S2/S3也已附原件和转录：S2百万级CD45数不是QC；S3用最后gene列。
说明与来源SHA见audit_data/evidence/GSE245552_S2_S3_readme.txt。
可单独离线重跑：python audit_data/audit_GSE245552_supplement_tables.py

Tabula Sapiens v2正常参照的另一个默认离线检查：
python audit_data/audit_tabula_sapiens_support.py
只重读所附52,298行肝/大肠必要metadata，不含42MB全图谱metadata或表达。
24总供者中只有TSP14同时有肝/大肠，严格髓系分别338/41；不能称24对。

追加GSE203608次要外部候选的默认离线检查：
python audit_data/audit_GSE203608_metadata.py
只读两份必要metadata列快照，保留原始文件全部行；未读表达。
临床in-vivo15人24组织59661细胞，4PTLM对09/27/44/49；四对髓系均>=100。
唯一四组织09的NC髓系只有11，不能计作>=50的四组织支持。
后增146988行文件含旧126014行完整重复，禁止直接相加；所有培养/共培养排除。
原图1B四对均MSS，但49原图StageIII与LiverMet标记并存，保留分期/逐例时序待核。
