"""Reaggregate preserved corpus metadata offline; no source-project dependency."""
from pathlib import Path
import csv,json,re

root=Path(__file__).resolve().parent
prov=root/'provenance'
metrics=json.loads((prov/'corpus_metrics.json').read_text(encoding='utf-8-sig'))
with (prov/'root_markdown_classification.csv').open(encoding='utf-8-sig',newline='') as f:
    rows=list(csv.DictReader(f))
terms=metrics['new_root_markdown']['deterministic_filename_heuristic']['terms']
total=sum(int(r['lines']) for r in rows)
assert (len(rows),total)==(315,77375)
variants={}
for label,excluded in [('original',[]),('without_MATRIX',['MATRIX']),
                       ('without_MATRIX_MAINLINE_CONSTITUTION',['MATRIX','MAINLINE','CONSTITUTION'])]:
    pattern='|'.join(t for t in terms if t not in excluded)
    selected=[r for r in rows if re.search(pattern,r['relative_path'],re.I)]
    lines=sum(int(r['lines']) for r in selected)
    variants[label]={'files':len(selected),'lines':lines,'file_percent':100*len(selected)/len(rows),
                     'line_percent':100*lines/total}
assert variants['original']['files']==162 and variants['original']['lines']==43603
assert variants['without_MATRIX']['files']==157 and variants['without_MATRIX']['lines']==42433
with (prov/'path_hash_relations.csv').open(encoding='utf-8',newline='') as f:
    relations=list(csv.DictReader(f))
counts={key:sum(r['relation']==key for r in relations) for key in ['identical','modified','missing','new']}
assert counts=={'identical':208,'modified':8,'missing':0,'new':580}
for r in relations:
    if r['relation']=='identical':assert r['baseline_hash']==r['continuation_hash']
    if r['relation']=='modified':assert r['baseline_hash']!=r['continuation_hash']
pdfs=json.loads((prov/'baseline_pdf_counts.json').read_text(encoding='utf-8'))
assert len(pdfs)==5 and sum(r['pages'] for r in pdfs)==85
result={'scope':'Preserved-metadata arithmetic only; historical PDF counting is inherited, not repeated here',
        'rows':len(rows),'lines':total,'path_relations':counts,'filename_sensitivity':variants,
        'baseline_manuscripts':5,'recorded_baseline_pages':85}
(root/'outputs').mkdir(exist_ok=True)
(root/'outputs/results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))
