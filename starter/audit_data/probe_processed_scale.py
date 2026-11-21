"""Numerically inspect a processed GEO matrix without inventing raw UMI counts.

Only 8 fixed columns are retained; the stream is read once across all gene rows.
This checks a candidate transformation, not a replacement for author provenance.
"""
from urllib.request import Request, urlopen
from pathlib import Path
import gzip, json, math, time, csv

URL='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164522/suppl/GSE164522_CRLM_PT_expression.csv.gz'
OUT=Path(__file__).resolve().parent/'geo'
start=time.monotonic(); k=8
rawsum=[0.0]*k; expsum=[0.0]*k; twosum=[0.0]*k; nnz=[0]*k
continuous=[0]*k; maxima=[0.0]*k; count=0; negative=0
try:
    with urlopen(Request(URL, headers={'User-Agent':'CRLM-data-audit/1.0'}),timeout=60) as response:
        with gzip.GzipFile(fileobj=response) as stream:
            header=stream.readline().decode().rstrip('\r\n')
            cells=next(csv.reader([header]))[1:k+1]
            for line in stream:
                prefix=line.decode().split(',',k+1)
                vals=[float(v.strip().strip('"')) for v in prefix[1:k+1]]
                count+=1
                for j,x in enumerate(vals):
                    rawsum[j]+=x; expsum[j]+=math.expm1(x)
                    twosum[j]+=2.0**x-1; nnz[j]+=(x!=0)
                    continuous[j]+=abs(x-round(x))>1e-6
                    maxima[j]=max(maxima[j],x); negative+=(x<0)
                if count%2000==0:
                    print(json.dumps({'rows':count,'elapsed_s':round(time.monotonic()-start,1)},ensure_ascii=False),flush=True)
    result={'status':'complete','url':URL,'n_genes_scanned':count,'n_columns_retained':k,
      'cell_columns':cells,'sum_as_published':rawsum,'sum_expm1':expsum,'sum_2power_minus1':twosum,
      'nonzero_features':nnz,'noninteger_features':continuous,'max_value':maxima,
      'negative_entries':negative,'elapsed_s':time.monotonic()-start,
      'warning':'Numerical compatibility is not verified author provenance; expm1 does not recover raw UMI counts.'}
except Exception as e:
    result={'status':'incomplete','url':URL,'rows_scanned':count,'elapsed_s':time.monotonic()-start,'error':str(e)}
(OUT/'GSE164522_processed_scale_probe.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
print(json.dumps(result,ensure_ascii=False),flush=True)
