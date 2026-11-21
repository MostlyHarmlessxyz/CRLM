# 文献与来源

完整文献表为 [bibliography.json](bibliography.json) 与 [reading_manifest.csv](reading_manifest.csv)。编号沿用原报告，准备优先级见 [论文与数据准备清单](../docs/PAPERS_AND_DATA_PREPARATION.md)。

[本地文件清单](local_files.csv) 记录项目中已归档的论文及附件。首批 12 篇核心论文已齐，包括 [M10（muscat）](papers/M10_2020.pdf)。B35 与 M24 共用 [scPAFA 正文](papers/M24_2024.pdf)，原有的 [B01 PDF](papers/B01_2020_classifier.pdf) 也列在清单中。

本地清单的 `archived_on` 为文件归档日期，`reference_url` 为文献入口。清单记录正文题名、DOI、版本、原文件名、页数、大小及 SHA-256。原研究的阅读深度见 `reading_manifest.csv`，起步包内已有的原表和审计证据见该表的 `starter_relative_paths`。

建议新文献记录：编号、标题、DOI/原始 URL、发布日期/版本、获取日期、文件相对路径和 SHA-256、全文/摘要/补充的阅读状态、关键页/表、直接用途和限制。补充表保留原名，放入 `supplements/<reference_id>/`；第三方代码记录 repository/commit/license/hash 放入 `code/<method>/`。同一 DOI 的两个研究用途编号不代表两篇独立文章。
