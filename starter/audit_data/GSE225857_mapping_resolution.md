# GSE225857 补充原表与患者代码核查

核查日期：2026-11-21。目标只限于核实论文六位患者与公开免疫元数据七个 `patients` 代码的关系；不据表达相似性、代码日期或排除后的样本数推断患者身份。

## 结论

**已从官方公开原始路径取得真实 Supplementary PDF 和 code.zip；补表提供了六位论文患者的临床、取样和治疗信息，但没有解释 `s0816` 的身份。患者映射问题仍未解决。**

Table S2、S3、S4 只列 `P0107`、`P0115`、`P0813`、`P0920`、`P1125`、`P1231` 六位 scRNA-seq 患者。Table S2 另列空间队列 `C1–C4`、`L1–L2`，不得当作这六位 scRNA-seq 患者的额外组织。没有任何一张表提到 `P0816` 或 `s0816`，也没有给出 `s0816 → s0813` 等别名表。

补充代码明确导入 `s0816` 数据，且以 barcode 前缀生成名为 `patients` 的变量。因此公开字段名称本身不足以证明有七位临床个体。不能把 `s0816` 合并进 `s0813`；不能称完整七代码数据为已确认的七患者队列。

## 1. 来源与完整性

原论文：Wang F et al. *Science Advances* 2023;9:eadf5464. DOI [10.1126/sciadv.adf5464](https://doi.org/10.1126/sciadv.adf5464)。

官方 BioStudies 记录：

- 记录 JSON：<https://www.ebi.ac.uk/biostudies/api/v1/studies/S-EPMC10275599>
- 公开文件目录元数据：<https://www.ebi.ac.uk/biostudies/api/v1/studies/S-EPMC10275599/info>，返回 `isPublic: true` 和官方 `httpLink`。
- PDF：<https://ftp.ebi.ac.uk/pub/databases/biostudies/S-EPMC/599/S-EPMC10275599/Files/sciadv.adf5464_sm.pdf>
- 代码：<https://ftp.ebi.ac.uk/pub/databases/biostudies/S-EPMC/599/S-EPMC10275599/Files/sciadv.adf5464_supplementary_code.zip>

独立官方备用入口也成功：<https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10275599/supplementaryFiles?includeInlineImage=false>。返回 ZIP 内的 PDF 和 code.zip 与 BioStudies 分别逐字节相同，核验记录为 `papers/GSE225857_mapping/cross_source_integrity.json`。

| 文件 | 字节数 | SHA-256 |
|---|---:|---|
| sciadv.adf5464_sm.pdf | 9,069,863 | `1250d41b389407f38e24f3ef1be7528d1375c1f61419a84b1abaf16607b9fd58` |
| sciadv.adf5464_supplementary_code.zip | 42,435 | `e683a84b93023247944fe260f29190c2b2f407cd7e22ba11d3a129aa723baf41` |

PDF 有 18 页，文件魔数 `%PDF-1.6`，能正常解析和渲染。临床表所在 PDF 第 13、14、15 页均已全文抽取，并逐页视觉核对。之前的 `adf5464_SM.failed.html` 仍是 HTML 失败页，不能充当附件。

## 2. 原表精确转录

完整转录保留原表行列与空白，位于：

- `papers/GSE225857_mapping/Table_S2.csv`：临床与取样，PDF 第 13 页；12 行，包括 6 位 scRNA-seq 和 6 位空间受试者。
- `papers/GSE225857_mapping/Table_S3.csv`：新辅助治疗，PDF 第 14 页；6 行。
- `papers/GSE225857_mapping/Table_S4.csv`：各部位 CD45 分选样本，PDF 第 15 页；6 行。
- `papers/GSE225857_mapping/tables_S2_S3_S4.json`：三表合并结构化记录、原表页码、来源和核验说明。

列标题只将换行合并，并将 CD45 上下标转成 `CD45+`/`CD45-`。空白没有改写为“无病变”或“缺少肿瘤”；原表说明仅为 `√` 表示存在且已取样。**三表没有 MSS/MSI/MMR 列，不能凭这些表补填该状态。**

| 原表患者 ID | 年龄 | 性别 | 原发部位 | S2 中取样部位 | S3 新辅助方案 | 周期 |
|---|---:|---|---|---|---|---:|
| P0107 | 75 | Male | Left colon | CN, CC, LN, LM, PB | FOLFOX | 3 |
| P0115 | 80 | Male | Left colon | CN, CC, LN, LM, PB | FOLFOX | 4 |
| P0813 | 58 | Female | Right colon | CN, CC, LN, LM, PB | FOLFOXIRI | 5 |
| P0920 | 64 | Male | Left colon | CC, LN, LM, PB | FOLFOX | 5 |
| P1125 | 41 | Female | Left colon | CN, CC, LN, LM | FOLFOXIRI+Bevacizumab | 8 |
| P1231 | 58 | Male | Left colon | CN, CC, LN, LM | FOLFOX | 3 |

原表六位 scRNA-seq 患者共计 27 个取样部位。S4 进一步明确：P1125 没有列出 CC/LM 的 CD45−细胞；P0920 没有 CN；P1125/P1231 没有 PB。公开 nonimmune 元数据恰含五个代码、十个代码×组织组，共 41,892 细胞；这与 S4 的样本模式一致，但不能据此解释免疫部分额外的 `s0816`。

## 3. 与公开免疫元数据的差异

下表来自本地已审计的 GEO 元数据汇总 `geo/GSM7058754_immune_meta_patient_tissue.csv`。公开代码与 P 代码数字后缀及取样模式对应，是可明确报告的文件间一致性；**原表未提供正式 crosswalk，额外代码的归属仍未知**。

| 公开 patients 代码 | 代码×组织组数 | 免疫细胞数 | 对应数字后缀是否在 S2/S3/S4 中 |
|---|---:|---:|---|
| s0107 | 5 | 34,980 | P0107 存在 |
| s0115 | 5 | 35,248 | P0115 存在 |
| s0813 | 5 | 27,597 | P0813 存在 |
| s0816 | 3 | 8,416 | **不存在** |
| s0920 | 4 | 26,229 | P0920 存在 |
| s1125 | 4 | 35,134 | P1125 存在 |
| s1231 | 4 | 28,869 | P1231 存在 |
| 合计 | 30 | 196,473 | 7 个代码，临床个体数不能由此直接确定 |

`s0816` 的三组为 LCL 5,334、LNL 615、PBL 2,467。去除这三组后为 27 组、188,057 个免疫细胞。**27 组吻合原文不足以证明应当删除 `s0816`，因为原文报告的 196,473 细胞恰好等于包含它的全体公开免疫细胞。**

## 4. 原始补充代码的证据

官方 code.zip 不含 `sampleinfo.csv`，只有脚本引用该本地文件。`supplementary_code/data_import_and_filter.R`：

- 第 51–57 行组装 immune_merge，包括 `s0813_3_4`、`s0813_5_7` 和 `s0816_2_4`。
- 第 59–60 行组装 tumor_merge，包括 `s0813_1_2` 和 `s0816_1`。
- 第 115 行读取作者计算机上的 `sampleinfo.csv`；第 116–130 行用其 date、sampletag、no、cell_type 建立组织标签。
- 第 158 行令 `patients` 等于细胞列名按下划线分隔所得前缀；没有显示从日期/批次代码到临床患者的另一个映射。
- 全 ZIP 文本检索 `0813`、`0816` 未发现二者合并、互换或别名声明。

代码中纳入过 `s0816` 的非免疫原始库，但最终公开 nonimmune 元数据没有此代码。这是另一个需要解释的处理阶段差异，不能据此替作者补写患者身份。

## 5. 对主方案的影响

1. **数据可获得性和临床治疗信息已经改善。** 不应继续写“Supplementary PDF 无法取得”；应改成“已取得 S2–S4，具体治疗可核，但额外代码身份未解”。
2. **完整队列仍不能作为已确认患者独立的主外部验证。** 主路线继续使用 GSE245552 发现、GSE164522 四组织确认、GSE178318 方向复现；GSE225857 保持附加、条件性的验证角色。
3. 可预先登记一项**排除身份未明代码的敏感性分析**：只纳入六个在表中有相同数字后缀的代码，按公开编码成组，明确分析的是 188,057 个免疫细胞和 27 代码×组织组；同时保留代码级而非已验证个体级的限定。它不能替代缺失 crosswalk，也不能凭单个敏感性结果将整个队列升级为严格确认集。
4. 若未来获得作者正式 crosswalk 或原始样本 manifest，应检查所有码的患者归属、时间点和重复取样，再决定同一患者内如何合并、保留重复测量或排除。得到映射前，不合并 `s0816/s0813`，不将 `s0816` 视作第七位已确认独立患者。
5. 表 S3 支持“六位已列临床患者接受过上述新辅助化疗”；不支持把跨队列差异归因于治疗。治疗种类、队列、平台仍高度混杂，也不能将 P0813 的 FOLFOXIRI 方案赋给 `s0816`。

## 6. 检索边界

此次只使用无需登录的公开官方入口，没有绕过验证码、付费墙或访问控制，没有联系作者。PMC OA API 与 article/bin 路径在当前尝试返回 404；这不表示附件不存在，BioStudies 与 EuropePMC REST 两个合法公开入口已取得字节一致的附件。附件取得成功与身份映射成功是两个不同结论。
