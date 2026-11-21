from pathlib import Path
import csv,json,re,hashlib,collections
import openpyxl
BASE=Path(__file__).resolve().parent; OUT=BASE.parent
rows=[]
FIELDS=['seed_set_id','proposed_axis','gene_symbol_original','source_ref','source_doi','source_url','source_location','source_definition','source_species','original_disease','cell_scope','evidence_kind','completeness','original_direction','original_weight','recommended_use','caveat']
def add(setid,axis,genes,ref,doi,url,loc,definition,disease,cell,kind,complete,direction,use,caveat):
 for gene,position in genes:
  rows.append(dict(zip(FIELDS,[setid,axis,gene,ref,doi,url,loc+position,definition,'Homo sapiens',disease,cell,kind,complete,direction,'',use,caveat])))
# Complete author-published marker blocks. No gene aliases, weights or membership invented.
s4=BASE/'NIHMS1944744-supplement-TABLE_S4.pdf.txt'
author_blocks=collections.defaultdict(list)
for page_no,page in enumerate(s4.read_text().split('\f'),1):
 for line in page.splitlines():
  bits=line.split()
  if len(bits)==2 and bits[0].startswith(('Mac_','Mono_')): author_blocks[bits[0]].append((bits[1],f'; PDF p{page_no}'))
for block,axis in [('Mac_SPP1','R_related_not_identical'),('Mac_CXCL10','C_related_not_identical')]:
 genes=author_blocks[block]
 assert len(genes)=={'Mac_SPP1':164,'Mac_CXCL10':94}[block]
 assert len({g for g,p in genes})==len(genes)
 add('B19_S4_'+block,axis,genes,'B19','10.1126/science.ade2292','https://ftp.ebi.ac.uk/pub/databases/biostudies/S-EPMC/760/S-EPMC10755760/Files/NIHMS1944744-supplement-TABLE_S4.pdf','Table S4; '+block,'该作者巨噬细胞状态在Table S4中的全部marker行；原始表两列cg/gene','Head and neck squamous cell carcinoma; cross-cancer replication in paper','Human tumor-associated macrophages','author_published_cluster_marker_block','complete_selected_author_block','marker_enriched_in_named_state','candidate_seed_requires_local_coverage_and_freezing','R/C为本项目拟议对应，非原作者同名程序；Mac_SPP1不等于纯脂质/基质；Mac_CXCL10不等于纯IFN-gamma；完整block不等于经跨CRLM验证签名。原始符号保持不变，未设权重。')
wbpath=next((BASE/'B06_supplement').glob('*.xlsx'));w=openpyxl.load_workbook(wbpath,read_only=True,data_only=True)
wburl='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9811165/supplementaryFiles?includeInlineImage=false'
want={'Extracellular matrix organization':'R_matrix_overlap','HDL-mediated lipid transport':'R_lipid_overlap','MHC class II antigen presentation':'A_MHCII_overlap','Antigen processing-Cross presentation':'A_crosspresentation_overlap'}
for rownum,row in enumerate(w['S5'].iter_rows(values_only=True),1):
 term=str(row[0]).strip()
 if term not in want:continue
 genes=str(row[3]).split(';');assert all(genes)
 add('B06_S5_'+want[term],want[term],[(g,'') for g in genes],'B06','10.1158/1078-0432.CCR-22-2041',wburl,f'Supplemental Table 5; S5 row {rownum}; Genes column D',term+'：肿瘤巨噬DEG与Reactome_2016条目相交的作者报告成员','Colorectal liver metastasis versus non-tumor liver','Human macrophages','author_reported_enrichment_overlap','complete_reported_overlap_not_complete_pathway','DEG overlap enriched pathway; no per-gene fold-change in S5','candidate_overlap_seed_not_full_pathway','不是完整Reactome通路，不是独立构建/验证的signature；不据富集断言功能增强；R/A共享溶酶体蛋白酶须检查。')
add('B06_text_EcoTyper_examples','R_related_not_identical',[(g,'') for g in ['SPP1','GPNMB','APOC1','APOE','TREM2','CTSB','LGALS3','FTL']],'B06','10.1158/1078-0432.CCR-22-2041','https://pmc.ncbi.nlm.nih.gov/articles/PMC9811165/','Results; paragraph beginning As an alternative approach; Supplementary Fig S2','作者正文对EcoTyper肿瘤巨噬状态所举marker','Colorectal liver metastasis','Human tumor-associated macrophages','verified_text_marker_examples','incomplete_examples','enriched_in_tumor_macrophage_state','interpretation_anchor_not_ready_to_score','原文列表以etc结尾，未补凑；不是EcoTyper完整signature；脂质/基质/吞噬等混合。')
for name,genes,axis in [('CXCL9_TAM',['CXCL9','IL4I1'],'C_related_not_identical'),('ISG15_TAM',['ISG15','CXCL10','CXCL11'],'C_related_not_identical')]:
 add('B20_text_'+name,axis,[(g,'') for g in genes],'B20','10.1158/2159-8290.CD-23-1300','https://pmc.ncbi.nlm.nih.gov/articles/PMC11294822/','Results first macrophage-clustering subsection; Fig 1D','作者正文描述的命名状态marker','Breast cancer and colorectal cancer; integrated four human studies','Human tumor-associated macrophages','verified_text_marker_examples','incomplete_examples','enriched_in_named_state','interpretation_anchor_not_ready_to_score','不能用IL4I1单基因等同IFN或呈递；该研究蛋白IL4I1组不能分辨CXCL9与ISG15转录组状态。')
add('B36_THP1_mechanistic_endpoints','C_experimental_endpoints_not_signature',[(g,'') for g in ['SPP1','CXCL9','CXCL10','STAT1']],'B36','10.1136/jitc-2025-013697','https://pmc.ncbi.nlm.nih.gov/articles/PMC12829364/','Results; Fig5-related paragraph; Supplementary Figs S9B-C and S13','SPP1-silenced THP1 macrophages的已核实分子/信号读出','THP1 cell line model; study also includes mouse tumors and human CRC tissue','Human cell-line-derived macrophage model','experimental_targets_or_readouts','not_a_gene_signature','SPP1 perturbation; CXCL9/10 and STAT1/p-STAT1 response','functional_readout_reference_only','THP1不是患者TAM；p-STAT1是蛋白磷酸化，不可用STAT1转录代替；鼠BMDM基因未转换后冒充人源signature。')
# Full background sets are preserved separately: these are published human-derived reused gene lists, not the small R signature.
bg=[]
s6=list(w['S6'].iter_rows(values_only=True))
for ci,name,expected,origin,doi in [(0,'Foam cell signature',26,'Human carotid atherosclerotic plaques; Fernandez2019','10.1038/s41591-019-0590-4'),(1,'Scar-associated macrophage signature',745,'Human liver cirrhosis; Ramachandran2019','10.1038/s41586-019-1631-3')]:
 sub=[]
 for rownum,vals in enumerate(s6,1):
  if rownum>=6 and vals[ci] is not None:
   sub.append({'set_id':'B06_S6_'+('Foam' if ci==0 else 'SAM'),'gene_symbol_original':str(vals[ci]),'source_ref':'B06','source_doi':'10.1158/1078-0432.CCR-22-2041','source_location':f'S6!{"AB"[ci]}{rownum}','source_definition':name,'original_study':origin,'original_study_doi':doi,'species':'Homo sapiens','completeness':'complete_column_as_reproduced_by_Sathe_not_rederived_from_original','recommended_use':'noncancer_background_or_specificity_check_only','caveat':'宽列表含housekeeping/核糖体/线粒体/性别等成分；原表标题把2019研究写为2020；不得直接等同冻结R。'})
 assert len(sub)==expected and len({z['gene_symbol_original'] for z in sub})==expected
 bg+=sub
with (OUT/'literature_program_seeds.csv').open('w',newline='',encoding='utf-8-sig') as f:
 z=csv.DictWriter(f,fieldnames=FIELDS);z.writeheader();z.writerows(rows)
with (BASE/'B06_S6_background_sets.csv').open('w',newline='',encoding='utf-8-sig') as f:
 z=csv.DictWriter(f,fieldnames=list(bg[0]));z.writeheader();z.writerows(bg)
# Complete macrophage rows from Table S5, retaining signed correlations as source facts, not score weights.
corrows=[]
for page_no,page in enumerate((BASE/'NIHMS1944744-supplement-TABLE_S5.pdf.txt').read_text().split('\f'),1):
 for line in page.splitlines():
  if not line.startswith('Macrophages'):continue
  # One PDF token merges the gene with its signed correlation; preserve source boundary rather than dropping it.
  line=line.replace('MIR4435-2HG-0.4171','MIR4435-2HG -0.4171')
  v=line.split();assert len(v)==7,(page_no,line)
  corrows.append(dict(cell_type=v[0],gene_symbol_original=v[1],cor=float(v[2]),original_rank=int(v[3]),S=float(v[4]),PCER=float(v[5]),FDR=float(v[6]),pdf_page=page_no,source_ref='B19',source_table='Table S5',source_doi='10.1126/science.ade2292',interpretation='correlation_with_CS_polarity_not_causal_effect_not_scoring_weight'))
assert len(corrows)==561
with (BASE/'B19_S5_macrophage_correlations.csv').open('w',newline='',encoding='utf-8-sig') as f:
 z=csv.DictWriter(f,fieldnames=list(corrows[0]));z.writeheader();z.writerows(corrows)
counts=collections.Counter(r['seed_set_id'] for r in rows)
manifest={'created_utc':'2026-11-21','main_rows':len(rows),'main_unique_genes':len({r['gene_symbol_original'] for r in rows}),'sets':dict(counts),'background_rows':len(bg),'B19_macrophage_correlation_rows':len(corrows),'status':'literature_seeds_NOT_project_frozen_signatures','file_integrity':[]}
for p in [wbpath,BASE/'NIHMS1944744-supplement-TABLE_S4.pdf.txt',BASE/'NIHMS1944744-supplement-TABLE_S5.pdf.txt',BASE/'NIHMS1944744-supplement-TABLE_S4.pdf',BASE/'NIHMS1944744-supplement-TABLE_S5.pdf',OUT/'literature_program_seeds.csv',BASE/'B06_S6_background_sets.csv',BASE/'B19_S5_macrophage_correlations.csv']:
 manifest['file_integrity'].append({'path':str(p.relative_to(OUT)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(BASE/'seed_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
print(json.dumps({k:v for k,v in manifest.items() if k!='file_integrity'},ensure_ascii=False,indent=2))
