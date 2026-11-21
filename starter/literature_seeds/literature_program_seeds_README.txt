文献程序种子说明（2026-11-21）
状态：仅为研究候选材料；不是本研究冻结签名，不得直接用于宣称已发现/已验证CRLM特异程序。

一、交付与来源
1. literature_program_seeds.csv：310条“集合—基因”记录、293个原始基因符号，共10组。每行保留来源、DOI、原始图表位置、适用物种/癌种/细胞、证据类型、完整性及使用边界。R/C/A是本项目拟议对应，并非作者原命名；不额外补基因、不指定新权重、不自动统一旧符号。
2. literature_program_sources/B06_S6_background_sets.csv：Foam 26 + scar-associated macrophage 745，完整照录Sathe补表S6两列，专用于非癌背景/特异性检查，不与主种子合并。
3. literature_program_sources/B19_S5_macrophage_correlations.csv：Bill补表S5中全部561条Macrophages记录，保留作者cor、rank、S、PCER、FDR与PDF页码。这些是CS极性关联统计量，不是本项目分数权重；不把正/负相关自动解释为C/R或抗/促肿瘤因果。
4. literature_program_sources/extract_seeds.py及seed_manifest.json：可复核提取脚本、计数、原始文件及CSV的SHA256。

二、哪些确有完整原始列表，哪些只是示例
B19 Bill 2023 Science，DOI 10.1126/science.ade2292：
- Table S4 Mac_SPP1整块164基因；Mac_CXCL10整块94基因。已取得原PDF，并抽看第9、27、30页核验表格与提取。
- 这两组是HNSCC作者定义的状态marker块；Mac_SPP1包括缺氧等成分，不能等同纯脂质/基质R；Mac_CXCL10不能等同纯IFN-gamma程序。它们也不是作者CXCL9:SPP1二基因极性比值的同义替代。
- Table S5 561个巨噬细胞基因的完整关联记录单独提供；其中CXCL9/10/11、STAT1/GBP及部分HLA/CD74均可直接核查，但没有将它们人为拼成一个新签名。
- 全文指出不同巨噬状态可以同时表达CXCL9与SPP1；Mac_APOE也更偏CXCL9侧，故不能仅依据“APOE=脂质”把所有APOE状态硬分到CS低侧。
- Table S4原PDF有其他非选中群中疑似Excel日期化符号（例如Mono_PADI4的Mar-01）。本次只纳入上述两块，其解析未发现日期形式；仍须在正式冻结前与目标数据基因标识核对，保留原文，不猜测修正。

B06 Sathe 2023 Clinical Cancer Research，DOI 10.1158/1078-0432.CCR-22-2041：
- 官方S1–S8 XLSX已获取，原文件在literature_program_sources/B06_supplement/。
- S5第23行“Extracellular matrix organization”：CD151, SDC2, FN1, P4HA1, CTSL, BSG, TIMP2, SPP1, ADAM9, ADAM8, P4HB, A2M, PPIB, CTSD, CTSB。
- S5第29行“HDL-mediated lipid transport”：ABCA1, APOE, A2M, PLTP。
- S5第18行“MHC class II antigen presentation”：CTSA, HLA-DMB, CTSL, CTSH, DYNLL1, CTSD, AP2M1, LGMN, CTSB。
- S5第36行“Antigen processing-Cross presentation”：PSMD8, PDIA3, UBB, CTSL, PSMB1, CALR, FCGR1A。
以上四项只是作者肿瘤巨噬DEGs与Reactome_2016条目的“完整报告交集”，不是完整Reactome通路，也不是独立训练/冻结的signature。S5不提供每一基因的方向/权重，因此CSV没有伪填fold-change或评分权重。需要方向应回到同一XLSX的S4并明确对照，而非根据通路名称推定。
- 正文EcoTyper段落的SPP1/GPNMB/APOC1/APOE/TREM2/CTSB/LGALS3/FTL仅是以etc.结束的8个例子，CSV明确不完整；未假装取得完整EcoTyper载荷。
- S6两列的原背景分别是人颈动脉粥样硬化（Fernandez 2019, PMID31591603, DOI10.1038/s41591-019-0590-4）与人肝硬化（Ramachandran 2019, PMID31597160, DOI10.1038/s41586-019-1631-3）。S6原标题把年份写为2020，正确发表年均为2019。此处保证Sathe所用列表的忠实照录，未声称重新从两篇最初数据推导。
- SAM745很宽，含核糖体、线粒体、免疫球蛋白、性别等成员，不能因为列名就当作精简巨噬纤维化签名；特别不能把GSE245552缺失的MT/RPL等基因补0。

B20 Matusiak 2024 Cancer Discovery，DOI10.1158/2159-8290.CD-23-1300：
- 正文Fig1D相关段明确CXCL9+TAM示例CXCL9/IL4I1，ISG15+TAM示例ISG15/CXCL10/CXCL11。它们是人乳腺癌/CRC四队列整合状态的少数marker，不是完整状态基因列表。
- 抗体定义的IL4I1组不能区分CXCL9与ISG15两类转录状态；IL4I1单基因也不等于有效抗原呈递或抗肿瘤功能。

B36 Wang 2026 JITC，DOI10.1136/jitc-2025-013697：
- SPP1/CXCL9/CXCL10/STAT1仅记录为人THP1来源巨噬细胞的干预/读出先例（Results Fig5相关段、Supp FigS9B-C/S13），不作为基因签名。
- 大量机制试验为鼠BMDM/鼠肿瘤；未把鼠基因转大写冒充人源列表。THP1不等于患者TAM，p-STAT1蛋白磷酸化也不能由STAT1转录量替代。

三、对第1周实施的具体用途
1. 导入后先按seed_set_id/evidence_kind筛选；functional_readout_reference_only和marker_examples不要自动送入主通路评分。完整原始列表也需先过测量覆盖、原始符号与当前符号映射审计。
2. 保存每队列每集合的可测交集、缺失基因数、版本、共同universe与缺失原因；GSE245552肿瘤特征过滤不能当作生物下调。不得按外部队列显著性删增成员。
3. R/A并非互斥：上述ECM与MHC-II交集共享CTSL/CTSD/CTSB，cross-presentation又含CTSL；应报告共享基因/溶酶体成员驱动情况及预设去共享成员敏感性。不能因名为“MHC-II呈递”便宣称患者实际呈递能力增强。
4. C的chemokine/IFN成分、A的HLA/加工/溶酶体成分和R的脂质/基质/缺氧成分需要分别检查。分离规则在发现集内完成并登记，再在确认集按冻结规则计算；这些文件只提供可追溯起点。
5. 如需最小基线，优先把作者完整Mac_CXCL10/Mac_SPP1块与明确的S5交集作为不同类型参照分别跑覆盖/稳定性检查，保留原始定义；不将全部集合并集后命名为文献验证R/C/A。

四、合法公开获取与核验
B19常规EuropePMC全文/补充接口当次返回500/非OA，未把错误XML当ZIP。随后经BioStudies官方公开记录（info明确isPublic=true）取得作者存档补充PDF：
https://www.ebi.ac.uk/biostudies/api/v1/studies/S-EPMC10755760
https://ftp.ebi.ac.uk/pub/databases/biostudies/S-EPMC/760/S-EPMC10755760/Files/NIHMS1944744-supplement-TABLE_S4.pdf
https://ftp.ebi.ac.uk/pub/databases/biostudies/S-EPMC/760/S-EPMC10755760/Files/NIHMS1944744-supplement-TABLE_S5.pdf
B19正文经NCBI官方BioC接口取得：
https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC10755760/unicode
B06补充下载：
https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9811165/supplementaryFiles?includeInlineImage=false
未突破登录/付费限制。原始PDF、XLSX与提取脚本均保留，失败响应单独保留供追溯。

五、起步包中的离线复核
本目录附两份原始事实型基因表PDF、其pdftotext输出和Sathe原始补充XLSX；未附论文全文或补充图片。literature_program_sources/extract_seeds.py只读这些本地小文件，额外需要openpyxl（本轮3.1.5），可重建三份CSV及哈希清单。正文示例部分是已核实原文的手工录入，脚本不会自动重新检索正文。原始符号与原始列表均保持，不做现代别名映射。
