CRLM 统计与数据起步包
版本：2026-11-21；Python 3.12.14 环境验证

用途
----
这是一套小型、可运行的患者级数据/统计起步实现，服务于“固定候选程序的配对效应、
四组织参照与组成/状态分解”。另附接收冻结细胞表示的透明线性attention MIL与mean-pooling
对照。它没有运行真实CRLM模型，不包含原始测序处理、细胞注释、UCell打分、负二项差异表达、
scGPT编码或完整tanh-ABMIL训练。synthetic_demo全部是合成数字。

安装与快速检查
--------------
在本目录运行（每行独立命令）：
    python -m pip install -r requirements.txt
    python selfcheck.py
    python mil_selfcheck.py
    python registered_tests_selfcheck.py

依赖固定为本轮实际验证版本；推荐Python 3.12。selfcheck_results.json记录环境与检查结果。
核心统计/MIL命令不需要GPU，也不联网。audit_data复核入口默认离线；显式选择--force的
完整矩阵复扫会重新下载公开数据，需另留带宽和时间；它不会替你申请受控原始测序数据。
mil_selfcheck_results.json记录独立的梯度/权重/角色隔离/单袋推理检查。
registered_tests_selfcheck_results.json记录冻结检验匹配、Holm、缺失家族成员与方向检查。

运行示例
--------
    python validate_manifest.py --manifest synthetic_demo/sample_manifest.csv --outdir demo_output/manifest
    python make_splits.py --donor-manifest demo_output/manifest/donor_manifest_validated.csv --n-folds 2 --output demo_output/folds.json
    python paired_effects.py --manifest demo_output/manifest/sample_manifest_validated.csv --scores synthetic_demo/program_scores.csv --outdir demo_output/effects
    python composition_decomposition.py --input synthetic_demo/composition.csv --output demo_output/composition.json
    python loss_counterexample.py --output demo_output/loss_counterexample.json
    python adjust_registered_tests.py --plan synthetic_demo/registered_tests.csv --summaries synthetic_demo/registered_effects_summary.csv --outdir demo_output/registered_tests

1. validate_manifest.py
-----------------------
输入CSV每行一个sample_uid。必需字段：
    dataset, donor_uid, sample_uid, tissue, lesion_uid, role,
    input_scale, patient_mapping_verified

tissue仅允许PT（原发肿瘤）、LM（肝转移灶）、NC（邻近结直肠）、NL（邻近肝）。
role仅允许development、external、exploratory、excluded。同一真实患者所有组织/病灶
必须共享同一个global donor_uid及同一个role。sample_uid、lesion_uid也要全局明确。
不同数据中都叫P1并不代表同一人；反之同一批生物样本再发表也不能因换accession当独立人。
跨来源身份审计后应将确定同人的样本映射到同一canonical donor_uid；仅加GSE前缀不足以防重。
同一患者多时间点/治疗前后不能自动视为同一配对窗口；本包假设输入已经按预定义临床窗口筛好。

input_scale允许raw_counts、normalized、log_normalized、corrected、unknown。
patient_mapping_verified只接受true/false或1/0，意为已由原始来源核实真实患者映射，
不是“ID字符串看起来能配对”。未核身份保留但eligible_for_paired_analysis=false。

可选字段：
    source_audit_scale: 来源审计记录的矩阵尺度。
    counts_provenance: 可追溯的原始counts来源/生成步骤说明。
    counts_values_verified: 是否已另外核查矩阵值域及来源（true/false）。
    can_nb_pseudobulk: 如果提供，脚本校验此声明；输出总会重新计算该能力。

只有input_scale=raw_counts、counts_provenance非空且counts_values_verified=true，
才输出can_nb_pseudobulk=true。已知processed来源与raw声明矛盾会报错；处理后表达声称
NB资格也报错。脚本不会打开表达矩阵，也无法证明人填写的来源是真实的；这不能代替数据审计。
整数dtype不证明原始counts，浮点dtype也不自动否认原始counts；不要sum(logexpr)、
expm1后四舍五入或把校正矩阵强行用于edgeR/DESeq2。

本轮GSE245552是一个具体反例：公开33份PT/LM计数矩阵各29569基因、6份NC/NL各29922，
前者是后者的严格子集，缺少353基因，其中有全部13个MT-基因和大量核糖体相关基因。
这批矩阵的整数性不等于完整未筛选CellRanger输出。肿瘤的空MT集合不能记成percent.mt=0，
不能因此判定高质量，也不能直接套用完整基因矩阵的UMI/gene阈值。PT/LM共同可测背景可用于
配对counts分析；四组织或跨队列编码/程序评分须冻结共同可测基因、记录覆盖并做敏感性分析。
未提供的基因不能补成生物学表达0。输入embedding的本包无法补回这些信息。

输出：sample_manifest_validated.csv、donor_manifest_validated.csv、manifest_validation.json。
重复sample、同患者跨role、同lesion跨患者/组织等硬错误以exit code 2退出。

2. make_splits.py
----------------
从上述donor表生成确定性的外层CV JSON；默认5折，demo为2折。必须2<=折数<=development
核实患者数。development中有未核身份则拒绝分折。external/exploratory/excluded均不进CV。
JSON逐折列出train/test donor和全部sample；四组织、多个病灶、技术重复始终跟随患者。
--seed控制打乱。留一患者可以设n-folds等于development患者数。

这里只产生外层分割；HVG/PCA/选基因/新程序拟合/超参数选择/早停仍需在训练折内完成。
用全部开发集先发现程序，再跑同一开发集CV，不会因有此JSON而自动免除筛选泄漏。
外部效应看过后更换程序或模型，需建立新探索性版本并撤回相应外部确认身份。

3. paired_effects.py
--------------------
scores CSV必需：sample_uid,lineage,program,score。可选n_cells仅作输入记录，永不作为患者权重。
其他样本信息统一来自manifest，不能在scores另放一套相互矛盾的donor/tissue列。
每个sample×lineage×program最多一行。score为空表示未可靠测量；0表示真实已测得的0。
输入应该是基因集合、方向、评分公式、亚型映射、覆盖阈值已经冻结后的样本程序读数。
本工具不会选择最佳候选，也不会把输入中细胞数少的亚型自动变成可靠表达估计。

顺序：同一病灶的技术样本取均值 -> 同患者同组织各病灶等权取均值 -> 每患者形成：
    D_T = LM - PT
    D_N = NL - NC
    I   = D_T - D_N
I只在同一患者四组织都有合格读数时计算；不会借别人的正常组织，不会缺失填0。
另输出D_T_Q、D_N_Q及complete_four_tissue：Q是同一批四组织合格患者。
背景图必须并列D_T_Q、D_N_Q、I；不可拿全部肿瘤配对的D_T均值减另一批正常配对的D_N均值。
将D_T与I比较时，另行报告同一四组织合格子集的D_T/D_N/I；不能用不同患者集合的均值相减。
每dataset、role、lineage、program独立汇总，所有患者等权，不把多个病灶变成多个N。
输入有未核身份/排除患者默认报错；如确需跳过，--skip-ineligible显式启用并写审计日志。

输出：lesion_means.csv、patient_tissue_means.csv、patient_contrasts.csv、effects_summary.csv、
analysis_audit.json。完全缺失但已提交的程序仍有N=0的汇总，便于看到失败/覆盖不足。
即使所有输入程序均不可测，也保留完整列结构，效应、原始p及区间为空；可直接接入登记工具。
均值差的患者bootstrap默认10000次，--seed固定；95%百分位CI在极小N下是近似精度，
不能称精确覆盖。N=1不输出bootstrap CI。细胞、病灶、重复抽袋都不单独重抽样为患者。

双侧sign test精确计算正向概率是否为0.5；符号CI使用Clopper–Pearson。
精确0差值为ties并从sign-test N中排除。全部ties给N=0、p=1、方向CI=[0,1]：表示无方向
信息，不是复现成功。ties半分率仅作描述，不将半个成功代入二项检验。
--tie-tolerance默认0；如因测量分辨率改容差，必须在看结果前规定并记在配置中。
输出正向比例指LM>PT（或所列对比>0）；预期负方向候选应按预先冻结方向解释。
这里输出未校正双侧p；多重比较族须在分析前定义，不自动跨所有组织/程序混作一个族。
随后将effects_summary.csv和已冻结的计划检验清单交给第8节的adjust_registered_tests.py。
I是观测的器官×肿瘤交互对比，不是转移因果效应。

4. composition_decomposition.py
-------------------------------
输入每行一个已核实患者×组织×亚型×固定程序：
    dataset, donor_uid, lineage, program, tissue, subtype, proportion, mean_score
本入口只比较PT/LM，proportion在每患者组织内需和为1。每种类型建议明确列出两组织行。
某组织亚型比例为0时，mean_score必须留空；将不存在的条件均值写0会报错。
如果比例>0但该类型均值不可靠，同样留空。这个脚本只核查代数输入，身份需先经manifest审计。

全共同支持时计算：
    D = sum(pi_LM*mu_LM - pi_PT*mu_PT)
    composition = sum((pi_LM-pi_PT)*(mu_LM+mu_PT)/2)
    state = sum((pi_LM+pi_PT)/2*(mu_LM-mu_PT))
它们满足D=composition+state，是线性均值恒等式，不是因果分解。
缺失一侧亚型均值时拒给完整分量，同时输出可测总体差（如可识别）、共同支持内重归一化
结果及两侧保留/缺失质量。共同支持结果的估计对象已改变，不能当作完整总体结果。
同一亚型在两组织都占0可排除，不影响任何总体均值。
这套公式不能直接分解注意力权重、非线性logit或sigmoid概率。

若以病灶等权为总体目标，应先形成一致的患者比例与条件均值：
    pi_patient,c = mean_lesion(pi_lesion,c)
    mu_patient,c = sum_lesion(pi_lesion,c*mu_lesion,c) / sum_lesion(pi_lesion,c)
分母为0则均值未知。不能分别对比例和条件均值都简单平均后，声称还等于病灶等权的总体均值。

5. loss_counterexample.py
-------------------------
固定肿瘤分数，只反转正常组织分数，展示原差分排序损失仍可获利。小平方惩罚不足以根治。
去掉差分奖励、单独使用对称正常差惩罚，没有该方向奖励，但也不保证消除器官混杂。
推荐主方案将背景检验外置；本脚本没有训练新方法或验证任何机制。

6. linear_attention_mil.py：冻结表示上的透明起步模型
--------------------------------------------------
每个输入NPZ对应一个固定谱系。每个key严格为sample_uid；每个value是该样本的
cells×features数值矩阵，可变细胞数、固定特征维度，不含object/pickle数组。
NPZ本身不带每细胞barcode。调用者须另外保存sample_uid、row_index、cell_uid的映射表及
SHA，确保导出的注意力数组与细胞行顺序对应；不能按另一次排序后的表达表直接解释权重。
训练先核验完整manifest，再仅读取选定development患者的PT/LM袋；NC/NL不进入BCE。
外部NPZ数组不会用于训练或标准化。缺失合格训练袋默认报错；确实不合格的袋需在manifest
事先声明embedding_eligible=false及非空embedding_exclusion_reason，不能看结果后静默丢人。

接口：
    python make_synthetic_embeddings.py
    python linear_attention_mil.py train --manifest synthetic_demo/sample_manifest.csv --embeddings synthetic_demo/cell_embeddings.npz --provenance synthetic_demo/embedding_provenance.json --lineage SyntheticImmune --pooling mean --standardize --output demo_output/mil_mean.npz
    python linear_attention_mil.py train --manifest synthetic_demo/sample_manifest.csv --embeddings synthetic_demo/cell_embeddings.npz --provenance synthetic_demo/embedding_provenance.json --lineage SyntheticImmune --pooling attention --standardize --output demo_output/mil_attention.npz
    python linear_attention_mil.py predict --model demo_output/mil_attention.npz --embeddings synthetic_demo/cell_embeddings.npz --sample E1_PT --output demo_output/mil_one_external.csv
    python linear_attention_mil.py predict --model demo_output/mil_attention.npz --embeddings synthetic_demo/cell_embeddings.npz --output demo_output/mil_all_predictions.csv --attention-output demo_output/mil_attention_weights.npz

在外层训练折拟合时再提供 --splits demo_output/folds.json --fold 0。训练器检查fold中患者
与当前development名单完全一致、train/test不交叉、每患者所有样本均被列入对应侧。
本模块不选择fold、超参数、候选或最优seed，不自动评估外部性能。

模型关系：
    mean:      z = mean_i(h_i), s = w^T z + b
    attention: a = softmax(H v), z = sum_i(a_i h_i), s = w^T z + b
它将标准ABMIL的v^T tanh(W h)注意力简化为h^T v，没有隐藏层W或tanh，不宣称复现完整ABMIL。
同维度D下，mean有D+1个参数，attention有2D+1个参数；训练输出记录参数数目。
注意力数组仅是归一化池化权重，没有细胞机制/因果/方向的自动解释。

训练目标（m=0）：
    L = sum_b(q_b * BCE(s_b,y_b)) + lambda_pair * mean_p softplus(-Delta_p)
        + lambda_l2/2 * (||w||^2 + ||v||^2)
其中mean模型没有v，b不受L2惩罚；q按患者→现有肿瘤部位→病灶→技术sample逐层等权。
Delta_p先同一病灶技术sample的logit等权，再同患者同组织病灶等权，最后LM减PT。
有BCE袋但无完整pair的患者只贡献BCE；lambda_pair>0却一个完整pair都没有时拒绝训练。
默认lambda_pair=1、lambda_l2=0.01只用于演示；真实值须在开发训练内规定/选择。
默认不拟合标准化。--standardize用上述q对训练袋的细胞均值/二阶差异加权计算中心和尺度，
保存在模型中，预测单袋时只应用保存参数，不利用其他待预测样本。HVG/PCA/scGPT等上游
步骤仍由调用者在正确训练折拟合或使用已审计的外部冻结权重；本模块无法挽救上游泄漏。

provenance JSON必需：
    lineage, encoder_id, representation_frozen, preprocessing_fit_scope,
    gene_universe_policy, gene_coverage_audited
两个布尔值必须是真正JSON true，不接受字符串。embedding_dimension可选但若提供须一致。
preprocessing_fit_scope允许external_pretrained、train_fold_only、cell_intrinsic、synthetic。
train_fold_only另需preprocessing_training_donors，且须等于当前预定训练患者名单。
external_pretrained必须意为整个任务表示/特征流程不使用外部任务数据拟合；若还有PCA/HVG
等任务适配，就不能只凭编码器被冻结而填这个值。建议额外记录checkpoint/词表/基因列表SHA、
code commit、gene coverage、预训练重叠状态、细胞纳入阈值与临床时间窗口。
gene_coverage_audited只是经过上游审计的声明；代码不读取原始基因，也不能证明声明真实。
GSE245552中正常与肿瘤的不同基因背景必须在编码前处理并记录，不能交给模型补零解决。

训练使用CPU float64与SciPy L-BFGS-B，固定seed。返回模型.npz、.training.json和.weights.csv。
不收敛时保存诊断JSON并报错，不写可预测模型。预测导出logit、probability_lm_label和cell数；
后者只是已获得组织的PT/LM标签概率，不是未来肝转移风险，正常组织上的概率也无临床解释。
was_training_sample仅检查sample ID是否相同，不是独立性证明。预测不接收真实标签。
实现将所需训练袋载入内存；它是小型可检查核心，不是百万细胞流式/GPU生产训练器。

验证：mil_selfcheck运行两种池化、两种pair系数的中心差分梯度检查，并验证病灶/患者权重、
完整manifest角色泄漏、错误fold、PCA训练范围、缺失袋、缺配对、单袋独立推理、cell顺序、
保存重载、固定seed，以及把外部数组改成NaN后仍能正常训练且训练参数不变。
合成embedding故意从标签生成信号，仅用于实现检查；任何演示拟合结果都不构成生物发现。

7. 方案冻结、候选登记与源码来源模板
---------------------------------
templates/analysis_plan.example.yaml是设计模板，不会自动改变CLI行为。它包含患者/病灶角色、
临床窗口、gene universe、平台覆盖、D_T与Q子集I、多重比较、缺失读数和scGPT逐细胞随机种子/
缓存身份。文件明确为synthetic_example=true、DRAFT_NOT_FROZEN；用于真实研究前需替换和冻结。

templates/candidate_registry.template.csv是空登记表；candidate_registry.synthetic.csv给出
与统计demo对应的Myeloid::fixed_demo、all_ties、missing_demo三个例子，其中后两者只作诊断，
不进入推断假设族。candidate_id用lineage::program；另锁定program_version。示例TOY_A/B/C
不是生物学基因，precomputed_synthetic_numbers不是UCell的真实运行输出。
candidate_gene_sets.synthetic.csv及frozen_gene_universe.synthetic.tsv的SHA写入登记表，
只是演示如何固定版本。coverage阈值是模板数值，不是本轮证明有效的推荐阈值。

measurement_coverage.template.csv用于记录每样本×候选的细胞/可测基因覆盖及缺失原因。
benchmark_run_registry.template.csv用于固定每次模型运行的输入、分割、覆盖、版本、耗时与失败记录。
external_access_log.template.csv用于披露每次外部访问和随后的方法修改；已知GSE164522的
作者标签/pilot访问已列为事件，不能将这种队列无条件称盲法确认。本轮已知访问另记录于
evidence/research_external_exposure_log.csv；只保留准确可核的日期及内容，未保留逐资源精确
时刻，不能补写虚构的时间。正式研究应从现在起逐次记录冻结版本和访问时间。

默认示例以固定候选×EXT测试列表构成PRIMARY_DT_EXTERNAL_V1并采用Holm；真实项目建议
主家族只含2—3个冻结候选×GSE164522，GSE178318预先列为次要方向/精度支持；背景I
另设冻结家族。各家族的错误率控制不等于全研究共同0.05控制，若需要后者须预先分配alpha
或定义统一家族。无法测量的计划测试保留在家族，原始p写缺失，若采用保守校正可在校正专用
向量中置p=1并明确标注，不能借缺失删掉难验证候选以缩小家族。

source_audit/source_provenance.json保存本轮相关官方源码入口、commit、文件SHA与核验范围：
  scGPT: cebd6fae655b9c585a4807daa3ac31bb764f06b4
  scMILD: eb1a284fd99b01fe340b0d6d1a934740012a9dae
  scMetaIntegrator: 97554c4c0eb7800864eb11a0b0af08f0c789895a
它们是读取时的版本，未声称复现原论文训练结果。scMILD默认先全量HVG再split，其sample ID
和单标签bag假设需要修改，NB自动编码器也不能直接吃164的processed表达。scGPT官方随机
ties/截断可能随batch和worker改变；缓存身份必须同时包含数据、基因范围、token规则、每细胞
seed和模型文件，而不能只以sample名或“scGPT”命名缓存。起步包不包含这些第三方训练代码。
完整方法文献及读取深度另见主报告附件methods_references.json；本包不宣称所有引用都全文精读。
已将scGPT 8个、scMILD 5个、scMetaIntegrator 4个所列文件与对应不可变commit重新下载的内容
逐一比较，SHA256全部一致。commit日期未知处保留null，不推测仓库版本就是论文实验版本。

8. adjust_registered_tests.py：严格匹配冻结检验与Holm校正
-------------------------------------------------------
计划CSV必需字段：
    test_id,family,cohort,role,lineage,program,contrast,expected_direction,
    direction_rule,alpha,test_method,adjustment,analysis_tier,freeze_status,
    freeze_timestamp_utc,freeze_version,synthetic_example

输入汇总CSV使用paired_effects.py的effects_summary.csv列名；cohort对应汇总中的dataset。
每个cohort×role×lineage×program×contrast只能对应一个计划检验，test_id也必须唯一。
计划或结果的重复键均报错，包括重复但未计划的结果。contrast支持D_T、D_N、D_T_Q、D_N_Q、I。
只接受test_method=exact_sign_two_sided和adjustment=Holm；还会由患者正/负方向计数核验
原始p确实为精确双侧sign test，因此必须保留p的完整数值精度，不能输入四舍五入的展示表。

expected_direction必须预设positive或negative。direction_rule当前仅支持
mean_and_sign_majority：同时要求患者平均差方向及非零差患者多数方向符合预期。
这是本起步包明确登记的保守方向规则；均值方向与多数方向也各自输出，便于发现不一致。
全ties的方向不可评价；p=1不意味着复现。N=0或整行结果缺失均保留raw_p为空，
只在adjustment_only_p列置1，仍占据原计划家族的一个名额。不能通过删除失败测量缩小家族。

同一个family的alpha、方法、层级、冻结版本和时间必须一致。每个家族独立进行Holm，
输出family_size_planned及其已测/缺失数。primary、secondary、background要求role=external；
开发结果仅可进入diagnostic家族。不同家族各用0.05并不等于全研究错误率共同0.05。
如希望全研究控制，应在访问结果前定义统一家族或分配alpha。

默认真实设计草稿templates/registered_tests.real_design_draft.csv有3个占位候选：
  PRIMARY_DT_164_V1：2—3个冻结候选×GSE164522的D_T，是主家族。
  SECONDARY_DT_178_V1：相同候选×GSE178318的D_T，提供次要方向/精度支持。
  BACKGROUND_I_164_V1：相同候选×GSE164522四组织Q子集的I，另作背景家族。
  SECONDARY_DT_203608_V1：可选的相同候选×GSE203608的D_T，独立次要方向支持。
应在发现阶段决定最终保留2个还是3个候选；不能读完外部结果后删掉候选或改变家族。
真实草稿以design_optional_before_freeze=true标记203608的三行；是否采用该队列必须在冻结
前决定。若不采用，应在冻结前整组移出计划；若保留，则其全部选定候选都计入独立次要家族，
不能看结果后删行。这个备注列不会让程序自动跳过行。未预先登记的203608结果只能单列探索。
203608只有4对PT–LM，不能把它加入164主家族、与178绝对分数混池或当作充分四组织确认。
178/203608的次要校正结果不自动成为主要复现的额外显著性门槛。背景需要同时展示同一Q患者的
D_T_Q、D_N_Q与I；单看I的p不足以区分所有机制，也不能证明“没有背景差异”。

草稿含REQUIRED占位、空冻结时间与DRAFT_NOT_FROZEN，程序会拒绝执行；需先完成真实设计。
templates/registered_tests.template.csv提供空表头。synthetic_demo/registered_tests.csv及
registered_effects_summary.csv是另外一组纯合成接口示例，含缺失、反向、全ties与计划外结果，
其SYN_GSE前缀不代表使用过真实队列。例子中的frozen及时间只演示字段，不是实际研究预注册。
脚本只能核查冻结声明的格式；真实先后关系须由外部访问日志、版本存档和SHA证据保证。

默认--unplanned-policy separate，将所有未计划结果写入unplanned_exploratory_results.csv；
它们永不进入冻结家族。若希望发现任何计划外结果即停止，可设置--unplanned-policy error。
输出另有registered_test_results.csv、registered_families.csv、registered_test_audit.json，
后者保存两份输入的SHA256。holm_reject与direction_matches分别输出；
statistical_direction_criterion_met只是二者同时成立，biological_replication_status始终是
not_assessed_by_this_tool。覆盖、病例独立、效应大小、组织背景和生物合理性仍须另外判断。

专用selfcheck当前32项全部使用合成数值，包括手算Holm、缺失维持家族大小、多个家族隔离、
重复拒绝、单侧p冒充双侧拒绝、相反方向、均值/多数方向分歧、全ties、全部缺失和实际CLI。
新增的整批不可测检查使用已有missing_demo，从paired_effects真实CLI连续流入登记CLI，
验证5个计划对比仍全部保留、raw p为空、校正专用p=1、家族大小仍为5。未访问真实外部效应。

验证覆盖与使用边界
------------------
selfcheck针对ID角色泄漏、冲突样本、raw/processed误声明、未核身份、病灶/患者等权、
缺失组织和完全缺失程序、全ties、精确小样本概率、组成共同支持及损失退化反例。
统计selfcheck实际调用五个CLI，MIL检查另实际调用train/predict，确保文件接口可衔接。
所有结果均来自合成数据。
本包支持实现起步与审计，不替代真实队列的数据来源核验、统计设计审查和生物学复现。

9. 原始文献程序种子与feature轴覆盖
--------------------------------
literature_seeds/literature_program_seeds.csv含310条集合×基因、293个原始符号，包括完整作者
marker块、富集交集、少数marker示例及实验读出；这些类型不能混为完整程序。B19的Mac_SPP1
与Mac_CXCL10来源于HNSCC，不是已经在CRLM验证的R/C；B06的S5列表是DEG与通路交集，
不是完整通路。先阅读同目录的literature_program_seeds_README.txt及逐行caveat。

    python literature_seeds/audit_seed_coverage.py

该命令只用标准库读取本包2份245代表性feature文件及4份335552面板，输出每原始符号覆盖、
每集合摘要及Jaccard共享基因表；没有表达、程序效应或外部标签拟合。缺少原符号可能是别名
或确实不在feature轴，不能全称生物下调。已存在的60行覆盖摘要也不是QC或生物验证结果。
若重新提取原事实型补表，可另安装openpyxl==3.1.5并运行：
    python literature_seeds/literature_program_sources/extract_seeds.py
该提取使用随附原始PDF的pdftotext输出及XLSX，保留原符号、顺序与SHA；不含论文全文，
marker示例是已核实正文的手工录入，不会在离线运行时再次检索来源。
