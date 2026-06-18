#!/usr/bin/env python3
"""
generate_report_json.py — Genomics report generator (sidecar JSON variant)

Variants are written as compact columnar JSON sidecars ({sid}.json),
fetched on demand via fetch(). Requires a local HTTP server to view:
    cd report/ && python -m http.server 8000

Usage:
    python generate_report_json.py <data_dir> [-o report/] [-r 0] [-c 1] [-f reference.fa] [--skip-json]
"""

import argparse, base64, gzip, json, mimetypes, re, sys
from datetime import datetime
from multiprocessing import Pool
from pathlib import Path

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(it, **kw):
        total = kw.get("total","?"); desc = kw.get("desc","")
        for i,item in enumerate(it,1):
            print(f"  {desc} [{i}/{total}]", file=sys.stderr); yield item

# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="Genomics report generator (JSON sidecar mode).")
    p.add_argument("data",                                help="Root data directory")
    p.add_argument("-o","--output-dir", default="report", help="Output directory (default: report/)")
    p.add_argument("-r","--max-vcf-rows", type=int, default=0,
                   help="Max VCF rows per sample, 0=unlimited (default: 0)")
    p.add_argument("-c","--ncpu", type=int, default=1,
                   help="Workers for parallel HTML rendering (default: 1)")
    p.add_argument("-g","--gff", default=None,
                   help="Sorted bgzipped+tabix-indexed GFF/GFF3 file (optional)")
    p.add_argument("-f", "--fasta", default=None,
                   help="Optional reference genome assembly FASTA (.fna/.fa)")
    p.add_argument("--skip-json", action="store_true",
                   help="Skip writing JSON sidecars (re-use existing files for fast UI testing)")
    return p.parse_args()

# ── Helpers ───────────────────────────────────────────────────────────────────

def _escape(s):
    return str(s).replace("&","&amp;").replace("<","&lt;").replace(">","&gt;")

IMAGE_EXTS = {".png",".jpg",".jpeg",".gif",".webp",".tiff",".bmp"}

def collect_images(subdir: Path) -> list:
    if not subdir.is_dir(): return []
    return sorted(p for p in subdir.iterdir() if p.suffix.lower() in IMAGE_EXTS)

def encode_image(path: Path) -> str:
    mime,_ = mimetypes.guess_type(str(path)); mime = mime or "image/png"
    return f"data:{mime};base64,{base64.b64encode(path.read_bytes()).decode()}"

def find_vcfs(sample_dir: Path) -> list:
    out = []
    for p in sorted(sample_dir.iterdir()):
        n = p.name
        if n.endswith(".vcf.gz") or n.endswith(".bcf.gz") or \
           n.endswith(".vcf")    or n.endswith(".bcf"):
            out.append(p)
    return out

# ── Markdown ──────────────────────────────────────────────────────────────────

def _md_inline(line: str) -> str:
    line = re.sub(r'`([^`]+)`', lambda m: f'<code>{_escape(m.group(1))}</code>', line)
    line = re.sub(r'\[([^\]]+)\]\(([^)]+)\)', r'<a href="\2">\1</a>', line)
    line = re.sub(r'\*\*(.+?)\*\*', r'<strong>\1</strong>', line)
    line = re.sub(r'__(.+?)__',     r'<strong>\1</strong>', line)
    line = re.sub(r'\*(.+?)\*',     r'<em>\1</em>',         line)
    line = re.sub(r'_(.+?)_',       r'<em>\1</em>',         line)
    line = re.sub(r'~~(.+?)~~',     r'<del>\1</del>',        line)
    return line

def parse_md(path: Path):
    if not path.exists(): return {}, ""
    text = path.read_text(errors="replace")
    fm, body = {}, text

    m = re.match(r"^---\s*\n(.*?)\n---\s*\n", text, re.DOTALL)
    if m:
        for line in m.group(1).splitlines():
            kv = line.split(":", 1)
            if len(kv) == 2: fm[kv[0].strip()] = kv[1].strip()
        body = text[m.end():]

    lines, html, i, in_ul, in_ol = body.splitlines(), [], 0, False, False

    def close_lists():
        nonlocal in_ul, in_ol
        if in_ul:  html.append("</ul>"); in_ul = False
        if in_ol:  html.append("</ol>"); in_ol = False

    while i < len(lines):
        raw = lines[i]
        line = raw.rstrip()
        if line.startswith("```"):
            close_lists()
            lang = line[3:].strip()
            code_lines = []
            i += 1
            while i < len(lines) and not lines[i].rstrip().startswith("```"):
                code_lines.append(_escape(lines[i]))
                i += 1
            lang_cls = f' class="language-{_escape(lang)}"' if lang else ""
            html.append(f'<pre><code{lang_cls}>{chr(10).join(code_lines)}</code></pre>')
            i += 1; continue
        if re.match(r'^[-*_]{3,}\s*$', line):
            close_lists(); html.append("<hr>"); i += 1; continue
        hm = re.match(r'^(#{1,6})\s+(.*)', line)
        if hm:
            close_lists()
            lvl = len(hm.group(1))
            html.append(f'<h{lvl}>{_md_inline(hm.group(2))}</h{lvl}>')
            i += 1; continue
        if line.startswith("> "):
            close_lists()
            html.append(f'<blockquote><p>{_md_inline(line[2:])}</p></blockquote>')
            i += 1; continue
        if "|" in line and i + 1 < len(lines) and re.match(r'^[\s|:-]+$', lines[i+1]):
            close_lists()
            headers = [c.strip() for c in line.strip().strip("|").split("|")]
            i += 2 
            html.append('<table class="md-table"><thead><tr>')
            for h in headers: html.append(f'<th>{_md_inline(h)}</th>')
            html.append('</tr></thead><tbody>')
            while i < len(lines) and "|" in lines[i]:
                cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                html.append('<tr>')
                for c in cells: html.append(f'<td>{_md_inline(c)}</td>')
                html.append('</tr>')
                i += 1
            html.append('</tbody></table>')
            continue
        if re.match(r'^[-*+] ', line):
            if in_ol: html.append("</ol>"); in_ol = False
            if not in_ul: html.append("<ul>"); in_ul = True
            html.append(f'<li>{_md_inline(line[2:])}</li>')
            i += 1; continue
        om = re.match(r'^(\d+)\. (.*)', line)
        if om:
            if in_ul: html.append("</ul>"); in_ul = False
            if not in_ol: html.append("<ol>"); in_ol = True
            html.append(f'<li>{_md_inline(om.group(2))}</li>')
            i += 1; continue
        if line == "":
            close_lists(); html.append(""); i += 1; continue
        close_lists()
        html.append(f'<p>{_md_inline(line)}</p>')
        i += 1

    close_lists()
    return fm, "\n".join(html)

# ── VCF parsing ───────────────────────────────────────────────────────────────

IMPACT_ORDER = {"HIGH":0,"MODERATE":1,"LOW":2,"MODIFIER":3,"":4}

def _parse_info(info_str):
    d = {}
    for part in info_str.split(";"):
        if "=" in part: k,v = part.split("=",1); d[k]=v
        else: d[part]=True
    return d

def _parse_csq_header(meta):
    for line in meta:
        m = re.search(r'##INFO=<ID=CSQ.*?Format: ([^"]+)"', line)
        if m: return [f.strip() for f in m.group(1).split("|")]
    return []

def _parse_ann_header(meta):
    for line in meta:
        m = re.search(r"##INFO=<ID=ANN.*?'([^']+)'", line)
        if m: return [f.strip().strip("'") for f in m.group(1).split("|")]
    return []

def _best(entries, key):
    return min(entries, key=lambda e: IMPACT_ORDER.get(e.get(key,""),4))

def _opener(p: Path):
    return gzip.open if p.suffix in (".gz",".bgz") else open

def read_vcf_header(vcf_path: Path):
    if not vcf_path.exists(): return [],[],[],[],[]
    meta,cols,samples = [],[],[]
    try:
        with _opener(vcf_path)(vcf_path,"rt",errors="replace") as fh:
            for raw in fh:
                line = raw.rstrip("\n")
                if line.startswith("##"): meta.append(line)
                elif line.startswith("#CHROM"):
                    cols = line.lstrip("#").split("\t")
                    # FIX: Truncate columns strictly to discard trailing duplicates
                    if len(cols) > 11:
                        cols = cols[:11]
                    if len(cols)>9:
                        samples = [s for s in cols[9:] if not s.startswith("2:")]
                    break
    except Exception as e:
        print(f"  [WARN] header read failed {vcf_path.name}: {e}", file=sys.stderr)
    return meta, cols, samples, _parse_csq_header(meta), _parse_ann_header(meta)

def stream_rows(vcf_path, col_headers, sample_names, csq_fmt, ann_fmt, max_rows=0):
    if not vcf_path.exists() or not col_headers: return
    raw_all = col_headers[9:]
    idx = {s: raw_all.index(s) for s in sample_names if s in raw_all}
    n = 0
    try:
        with _opener(vcf_path)(vcf_path,"rt",errors="replace") as fh:
            for raw in fh:
                line = raw.rstrip("\n")
                if line.startswith("#"): continue
                parts = line.split("\t")
                if len(parts)<8: continue
                chrom,pos,_,ref,alt,qual,filt,info_str = parts[:8]
                scols = parts[9:] if len(parts)>9 else []
                
                info = _parse_info(info_str)
                info_dp = info.get("DP", "0") # Extract global Depth
                
                format_col = parts[8].split(":") if len(parts) > 8 else []
                dp_idx = format_col.index("DP") if "DP" in format_col else -1

                gt = {}
                for s,ci in idx.items():
                    if scols and ci < len(scols):
                        sdata = scols[ci].split(":")
                        call = sdata[0] if sdata else "./."
                        depth = sdata[dp_idx] if dp_idx != -1 and dp_idx < len(sdata) else "?"
                        gt[s] = f"{call} ({depth}x)"
                    else:
                        gt[s] = "./. (?x)"
                
                vep, ann = [], []
                if csq_fmt and "CSQ" in info:
                    for blk in info["CSQ"].split(","):
                        vs=blk.split("|"); vep.append({f:(vs[j] if j<len(vs) else "") for j,f in enumerate(csq_fmt)})
                if ann_fmt and "ANN" in info:
                    for blk in info["ANN"].split(","):
                        vs=blk.split("|"); ann.append({f:(vs[j] if j<len(vs) else "") for j,f in enumerate(ann_fmt)})
                
                vi = _best(vep,"IMPACT").get("IMPACT","") if vep else ""
                ai = _best(ann,"Annotation_Impact").get("Annotation_Impact","") if ann else ""
                impact = min(vi,ai,key=lambda x:IMPACT_ORDER.get(x,4)) if (vi or ai) else ""
                
                yield {"chrom":chrom,"pos":pos,"ref":ref,"alt":alt,"qual":qual,"filter":filt,
                       "info_dp": info_dp, "vep":vep,"ann":ann,"lof":info.get("LOF",""),
                       "nmd":info.get("NMD",""),"impact":impact,"samples":gt}
                n+=1
                if max_rows and n>=max_rows: return
    except Exception as e:
        print(f"  [WARN] parse error {vcf_path.name}: {e}", file=sys.stderr)


def write_variant_json(vcf_path, json_path, col_headers, sample_names, csq_fmt, ann_fmt, max_rows=0):
    SCALAR = ["chrom","pos","ref","alt","qual","filter","info_dp","impact","lof","nmd",
              "vep_symbol","vep_consequence","vep_hgvsc","vep_hgvsp",
              "vep_aa","vep_biotype","vep_class",
              "ann_gene","ann_annotation","ann_hgvsc","ann_hgvsp","ann_biotype",
              "vep_n","ann_n"]
    all_fields = SCALAR + [f"s:{s}" for s in sample_names]
    n = 0
    coord_ranges = {} 
    with open(json_path,"w",encoding="utf-8") as out:
        out.write('{"fields":'); out.write(json.dumps(all_fields,separators=(',',':'))); out.write(',"rows":[')
        for v in stream_rows(vcf_path,col_headers,sample_names,csq_fmt,ann_fmt,max_rows):
            vb = _best(v["vep"],"IMPACT") if v["vep"] else {}
            ab = _best(v["ann"],"Annotation_Impact") if v["ann"] else {}
            row = [v["chrom"],v["pos"],v["ref"],v["alt"],v["qual"],v["filter"],v.get("info_dp","0"),v["impact"],
                   bool(v["lof"]),bool(v["nmd"]),
                   vb.get("SYMBOL",""),vb.get("Consequence",""),vb.get("HGVSc",""),
                   vb.get("HGVSp",""),vb.get("Amino_acids",""),vb.get("BIOTYPE",""),vb.get("VARIANT_CLASS",""),
                   ab.get("Gene_Name",""),ab.get("Annotation",""),ab.get("HGVS.c",""),
                   ab.get("HGVS.p",""),ab.get("Transcript_BioType",""),
                   len(v["vep"]),len(v["ann"]),
                   *[v["samples"].get(s,"./.")for s in sample_names]]
            if n>0: out.write(",")
            out.write(json.dumps(row,separators=(',',':')))
            n += 1
            try:
                pos, c = int(v["pos"]), v["chrom"]
                if c not in coord_ranges: coord_ranges[c] = [pos, pos]
                else:
                    if pos < coord_ranges[c][0]: coord_ranges[c][0] = pos
                    if pos > coord_ranges[c][1]: coord_ranges[c][1] = pos
            except: pass
        out.write("]}")
    return n, {c: (r[0], r[1]) for c, r in coord_ranges.items()}

# ── Discovery ─────────────────────────────────────────────────────────────────

MARKERS = {"genome_overview","nucleotide_view","sample.md"}

def _is_sample_dir(d: Path):
    try:
        names = {p.name for p in d.iterdir()}
        return bool(names & MARKERS) or any(".vcf" in n for n in names)
    except PermissionError: return False

def discover_samples(data_dir: Path):
    if not data_dir.is_dir(): return []
    if _is_sample_dir(data_dir): return [data_dir]
    out = []
    for d in sorted(data_dir.iterdir()):
        if (d.is_dir() or (d.is_symlink() and d.resolve().is_dir())) and _is_sample_dir(d):
            out.append(d)
    return out

# ── IGV Track Builder ─────────────────────────────────────────────────────────

def get_bam_locus(bam_path: Path) -> str:
    import subprocess
    try:
        idx_res = subprocess.run(["samtools", "idxstats", str(bam_path)], capture_output=True, text=True)
        best_chrom, max_reads = "", 0
        for line in idx_res.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) >= 3 and int(parts[2]) > max_reads:
                max_reads, best_chrom = int(parts[2]), parts[0]
        if not best_chrom: return "all"
        view_res = subprocess.run(["samtools", "view", str(bam_path), best_chrom], capture_output=True, text=True)
        first_pos = int(view_res.stdout.splitlines()[0].split("\t")[3])
        return f"{best_chrom}:{max(1, first_pos - 2000)}-{first_pos + 2000}"
    except Exception: return "all"

# ── CSS / JS (shared) ─────────────────────────────────────────────────────────

CSS = """
:root{--bg:#0f1117;--surf:#181c27;--surf2:#1e2235;--brd:#2a2f45;--txt:#d0d6f0;
  --muted:#6b7494;--acc:#4f8ef7;--acc2:#7c5be8;--r:6px;
  --mono:'JetBrains Mono','Fira Mono',monospace;--sans:'IBM Plex Sans','Segoe UI',system-ui,sans-serif;}
*{box-sizing:border-box;margin:0;padding:0}html{scroll-behavior:smooth}
body{background:var(--bg);color:var(--txt);font-family:var(--sans);font-size:14px;
  line-height:1.6;display:flex;min-height:100vh}
a{color:var(--acc);text-decoration:none}a:hover{text-decoration:underline}
/* sidebar */
#sidebar{width:220px;min-width:220px;background:var(--surf);border-right:1px solid var(--brd);
  position:fixed;top:0;left:0;bottom:0;overflow-y:auto;z-index:100;padding-bottom:2rem}
.sb-head{padding:1rem;border-bottom:1px solid var(--brd);background:var(--surf2);position:sticky;top:0;z-index:1}
.sb-head .h1{font-size:13px;font-weight:700;letter-spacing:.06em;text-transform:uppercase;color:var(--acc)}
.sb-head .sub{font-size:11px;color:var(--muted);margin-top:2px}
.sb-search{padding:.5rem 1rem;border-bottom:1px solid var(--brd)}
.sb-search input{width:100%;background:var(--bg);border:1px solid var(--brd);color:var(--txt);
  padding:.3rem .5rem;border-radius:var(--r);font-size:12px}
.sb-search input:focus{outline:1px solid var(--acc)}
.sb-idx{display:block;padding:.6rem 1rem;font-size:12px;color:var(--acc);
  border-bottom:1px solid var(--brd);font-weight:600}
.sb-idx:hover{background:var(--surf2)}
.sb-item{display:block;padding:.45rem 1rem;font-size:12px;color:var(--txt);
  border-left:3px solid transparent;transition:all .15s;
  white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
.sb-item:hover,.sb-item.active{background:var(--surf2);border-left-color:var(--acc);color:var(--acc);text-decoration:none}
.sb-item.current{border-left-color:var(--acc2);color:var(--acc2);background:var(--surf2)}
/* main */
#main{margin-left:220px;padding:2rem;max-width:1400px;width:100%}
/* index */
.idx-section h2{font-size:1.4rem;margin-bottom:1rem;border-bottom:1px solid var(--brd);padding-bottom:.5rem}
.idx-table{width:100%;border-collapse:collapse;font-size:13px}
.idx-table th{background:var(--surf2);padding:.5rem .75rem;text-align:left;
  border-bottom:1px solid var(--brd);color:var(--muted);text-transform:uppercase;font-size:11px;letter-spacing:.05em}
.idx-table td{padding:.5rem .75rem;border-bottom:1px solid var(--brd)}
.idx-table tr:hover td{background:var(--surf2)}
.tc{text-align:center}
/* sample page */
.back-link{font-size:12px;color:var(--muted);display:inline-block;margin-bottom:.5rem}
.back-link:hover{color:var(--acc)}
.sample-name{font-size:1.6rem;font-weight:700;font-family:var(--mono);margin-bottom:.5rem}
.meta-dl{display:grid;grid-template-columns:auto 1fr;gap:.2rem .75rem;font-size:12px;
  margin:.5rem 0;background:var(--surf2);padding:.5rem .75rem;border-radius:var(--r);
  border:1px solid var(--brd);max-width:600px}
.meta-dl dt{color:var(--muted);font-weight:600}.meta-dl dd{color:var(--txt)}
.sample-notes{margin-top:.75rem;background:var(--surf2);border-left:3px solid var(--acc2);
  padding:.75rem 1rem;border-radius:0 var(--r) var(--r) 0;font-size:13px;max-width:800px}
/* master toggle */
.master-bar{display:flex;gap:.5rem;align-items:center;margin:1rem 0;flex-wrap:wrap}
.master-btn{background:var(--surf2);border:1px solid var(--brd);color:var(--muted);
  padding:.3rem .75rem;border-radius:var(--r);font-size:12px;cursor:pointer;font-family:var(--sans)}
.master-btn:hover{border-color:var(--acc);color:var(--acc)}
/* collapsible sections */
.section{margin-bottom:1.5rem;border:1px solid var(--brd);border-radius:var(--r);overflow:hidden}
.section-hdr{display:flex;align-items:center;gap:.5rem;padding:.6rem .75rem;
  background:var(--surf2);cursor:pointer;user-select:none}
.section-hdr:hover{background:#232740}
.section-hdr h3{font-size:.95rem;font-weight:700;flex:1}
.section-arrow{font-size:.75rem;color:var(--muted);transition:transform .2s}
.section-arrow.open{transform:rotate(90deg)}
.section-body{display:none;padding:.75rem}
.section-body.open{display:block}
.count{color:var(--muted);font-size:.85em;font-weight:400}
.muted{color:var(--muted);font-style:italic;font-size:13px}


/* show help icons */
.th-wrap { 
  display: flex; 
  align-items: center; 
  gap: 6px; 
  position: relative; 
}
.tt-icon {
  display: inline-flex; align-items: center; justify-content: center;
  width: 14px; height: 14px; border-radius: 50%; background: #34495e;
  color: #ffffff; font-size: 10px; cursor: help; position: relative;
  font-family: sans-serif; font-weight: bold;
}
/* Drop tooltips downwards so they never clip out of frame */
.tt-icon::after {
  content: attr(data-tt); position: absolute; 
  top: 140%; left: 50%; transform: translateX(-50%);
  background: #2c3e50; color: #ffffff; padding: 8px 12px; 
  border: 1px solid #3498db; border-radius: 4px; font-size: 11px; 
  white-space: normal; width: max-content; max-width: 160px;
  z-index: 9999; box-shadow: 0 4px 12px rgba(0,0,0,0.3); 
  text-transform: none; font-weight: normal; opacity: 0; 
  pointer-events: none; transition: opacity 0.15s;
}
.tt-icon:hover::after { opacity: 1; }

/* gallery */
.gallery{display:flex;flex-wrap:wrap;gap:.75rem}
.img-figure{background:var(--surf2);border:1px solid var(--brd);border-radius:var(--r);
  overflow:hidden;cursor:zoom-in;transition:border-color .15s}
.img-figure:hover{border-color:var(--acc)}
.img-figure img{display:block;max-width:min(100%,900px);max-height:420px;object-fit:contain;background:#000}
.img-figure figcaption{padding:.3rem .5rem;font-size:11px;color:var(--muted);font-family:var(--mono)}
/* lightbox */
#lightbox{display:none;position:fixed;inset:0;background:rgba(0,0,0,.92);z-index:1000;
  justify-content:center;align-items:center;cursor:zoom-out;flex-direction:column;gap:.5rem}
#lightbox.open{display:flex}
#lightbox img{max-width:96vw;max-height:90vh;object-fit:contain;border-radius:var(--r)}
#lb-cap{color:#aaa;font-size:12px;font-family:var(--mono)}
#lb-close{position:fixed;top:1rem;right:1rem;background:rgba(255,255,255,.1);border:none;color:#fff;
  font-size:1.5rem;cursor:pointer;border-radius:50%;width:36px;height:36px;display:flex;align-items:center;justify-content:center}
#lb-close:hover{background:rgba(255,255,255,.2)}
/* load btn */
.load-btn{background:var(--surf2);border:1px solid var(--acc);color:var(--acc);
  padding:.5rem 1.25rem;border-radius:var(--r);font-size:13px;cursor:pointer;
  font-family:var(--sans);transition:background .15s}
.load-btn:hover{background:rgba(79,142,247,.15)}
/* filter bar */
.filter-bar{display:flex;flex-wrap:wrap;gap:.8rem;align-items:center;margin-bottom:.75rem;
  padding:.5rem .75rem;background:var(--surf2);border-radius:var(--r);border:1px solid var(--brd)}
.filter-bar span{color:var(--muted);font-size:12px;margin-right:.25rem}
.f-btn{background:var(--bg);border:1px solid var(--brd);color:var(--muted);
  padding:.2rem .6rem;border-radius:var(--r);font-size:11px;cursor:pointer;font-family:var(--sans)}
.f-btn.active{border-color:var(--acc);color:var(--acc);background:rgba(79,142,247,.12)}
.search-box{background:var(--bg);border:1px solid var(--brd);color:var(--txt);
  padding:.25rem .6rem;border-radius:var(--r);font-size:12px;width:220px;margin-left:auto}
.search-box:focus{outline:1px solid var(--acc)}
/* variant table */
.tbl-scroll{overflow-x:auto; max-width: 100vw;}
.var-table{width: max-content;border-collapse:collapse;font-size:12px;min-width:100%;table-layout: auto;}
.var-table th{background:var(--surf2);padding:10px 12px;text-align:left;
  border-bottom:1px solid var(--brd);color:var(--muted);font-size:11px;
  text-transform:uppercase;letter-spacing:.04em;white-space:nowrap;position:sticky;top:0;z-index:10;cursor:pointer;}
.var-table th:hover{color:var(--txt);}
.var-table td{padding:.4rem .6rem;border-bottom:1px solid var(--brd);vertical-align:middle; white-space:nowrap; max-width:250px; overflow:hidden; text-overflow:ellipsis;}
.var-row{cursor:default; transition: opacity 0.3s, filter 0.3s;}
.var-row:hover td{background:var(--surf2)}
.var-row[data-impact=HIGH]     td:first-child{border-left:3px solid #f44}
.var-row[data-impact=MODERATE] td:first-child{border-left:3px solid #f90}
.var-row[data-impact=LOW]      td:first-child{border-left:3px solid #48f}

/* New Tiers */
.var-row.tier-1 { opacity: 1; }
.var-row.tier-2 { opacity: 0.55; }
.var-row.tier-3 { opacity: 0.20; filter: grayscale(100%); }
.var-row.saved { background-color: rgba(46, 204, 113, 0.15) !important; opacity: 1 !important; filter: none !important; border-left: 3px solid #2ecc71 !important;}
.var-row.discarded { opacity: 0.10 !important; filter: grayscale(100%) !important; text-decoration: line-through; border-left: 3px solid #e74c3c !important;}

.sc{max-width:200px;white-space:normal;word-break:break-word;font-size:11px}
.gt{font-family:var(--mono);font-size:11px;text-align:center;color:var(--muted)}
.badge{display:inline-block;padding:.1rem .4rem;border-radius:3px;font-size:10px;
  font-weight:700;letter-spacing:.04em;margin-right:.25rem;font-family:var(--mono)}
.lof-tag,.nmd-tag{display:inline-block;background:#7c5be8;color:#fff;padding:.1rem .35rem;
  border-radius:3px;font-size:9px;font-weight:700;margin-right:.2rem;vertical-align:middle}
.nmd-tag{background:#2a9d8f}
code{font-family:var(--mono);font-size:12px;background:var(--surf2);padding:.1rem .3rem;border-radius:3px}
.viewer-link{display:inline-block;margin-left:.6rem;font-size:11px;color:var(--acc);border:1px solid var(--brd);padding:.2rem .5rem;border-radius:var(--r);text-decoration:none;vertical-align:middle}
.viewer-link:hover{border-color:var(--acc);background:rgba(79,142,247,.08)}
.ts-footer{font-size:10px;color:var(--muted);padding:.4rem 2rem 1rem;opacity:.6}
.md-table{width:100%;border-collapse:collapse;font-size:13px;margin:.75rem 0}
.md-table th{background:var(--surf2);padding:.4rem .75rem;text-align:left;border-bottom:2px solid var(--brd);color:var(--muted);font-size:11px;text-transform:uppercase;letter-spacing:.04em}
.md-table td{padding:.35rem .75rem;border-bottom:1px solid var(--brd)}
.md-table tr:hover td{background:var(--surf2)}
blockquote{border-left:3px solid var(--acc2);margin:.5rem 0;padding:.4rem .75rem;background:var(--surf2);border-radius:0 var(--r) var(--r) 0;font-size:13px;color:var(--muted)}
pre{background:var(--surf2);border:1px solid var(--brd);border-radius:var(--r);padding:.75rem 1rem;overflow-x:auto;margin:.5rem 0}
pre code{background:none;padding:0;font-size:12px}
hr{border:none;border-top:1px solid var(--brd);margin:1rem 0}
.variants-section{margin-top:1.5rem}
.variants-title{font-size:1rem;font-weight:700;margin-bottom:1rem;color:var(--txt)}
.vcf-block{margin-bottom:2rem;padding-top:.75rem;border-top:1px solid var(--brd)}
.vcf-name{font-size:.85rem;font-weight:700;font-family:var(--mono);color:var(--muted);
  margin-bottom:.5rem;letter-spacing:.03em}
/* GFF track */
.gff-sticky{position:sticky;top:0;z-index:40;background:var(--bg);
  border-bottom:1px solid var(--brd);padding:.4rem 0 .3rem;margin-bottom:.75rem}
.gff-header{display:flex;align-items:center;gap:.75rem;padding:0 .25rem .3rem;font-size:11px}
.gff-title{font-weight:700;color:var(--muted);text-transform:uppercase;letter-spacing:.05em}
.gff-coords{color:var(--muted);font-family:var(--mono);font-size:10px}
.gff-svg{display:block;overflow:visible}
#kbd-hint{position:fixed;bottom:.75rem;right:.75rem;font-size:10px;color:var(--muted);
  background:var(--surf);border:1px solid var(--brd);padding:.3rem .6rem;border-radius:var(--r);opacity:.7}
"""

JS = """
// ── global state ──────────────────────────────────────────────────────────
const _varData  = {};   // jsid → parsed data object
const _posIndex = {};   // jsid → {chrom: [[pos, rowIdx], ...]} sorted by pos
const _viewMode = {};   // jsid → 'window' | 'all'
const _currentIndices = {}; // jsid -> array of currently rendered rows
const _sortState = {};

// Global memory caches to preserve sort and triage states between redraws
window.sortStates = {};
window.triageRegistry = {};

function setTriageState(jsid, rowIdx, action) {
    const rowEl = document.getElementById(`vrow-${jsid}-${rowIdx}`);
    const key = `${jsid}_${rowIdx}`;
    
    if (!window.triageRegistry[jsid]) window.triageRegistry[jsid] = {};

    if (window.triageRegistry[jsid][rowIdx] === action) {
        // Toggle off if clicking the same action button again
        delete window.triageRegistry[jsid][rowIdx];
        if (rowEl) rowEl.className = 'var-row';
    } else {
        window.triageRegistry[jsid][rowIdx] = action;
        if (rowEl) {
            rowEl.className = `var-row ${action === 'save' ? 'saved' : 'discarded'}`;
        }
    }
}

// ── lightbox ──────────────────────────────────────────────────────────────
function openLightbox(fig){
  document.getElementById('lb-img').src = fig.querySelector('img').src;
  document.getElementById('lb-cap').textContent = (fig.querySelector('figcaption')||{textContent:''}).textContent;
  document.getElementById('lightbox').classList.add('open');
}
document.getElementById('lightbox').addEventListener('click',()=>document.getElementById('lightbox').classList.remove('open'));
document.getElementById('lb-close').addEventListener('click',e=>{e.stopPropagation();document.getElementById('lightbox').classList.remove('open');});
document.addEventListener('keydown',e=>{ if(e.key==='Escape') document.getElementById('lightbox').classList.remove('open'); });

// ── collapsible sections ──────────────────────────────────────────────────
function toggleSection(id){
  const body=document.getElementById('sb-'+id);
  const arrow=document.getElementById('ar-'+id);
  if(!body) return;
  const open=body.classList.toggle('open');
  if(arrow) arrow.classList.toggle('open',open);
}
function expandAll(){
  document.querySelectorAll('.section-body').forEach(b=>b.classList.add('open'));
  document.querySelectorAll('.section-arrow').forEach(a=>a.classList.add('open'));
}
function collapseAll(){
  document.querySelectorAll('.section-body').forEach(b=>b.classList.remove('open'));
  document.querySelectorAll('.section-arrow').forEach(a=>a.classList.remove('open'));
}

// ── sidebar filter ────────────────────────────────────────────────────────
function sbFilter(inp){
  const q=inp.value.toLowerCase();
  document.querySelectorAll('.sb-item').forEach(el=>{
    el.style.display=el.textContent.toLowerCase().includes(q)?'':'none';
  });
}

// ── badge / esc helpers ───────────────────────────────────────────────────
const IC={HIGH:['#f44','#fff'],MODERATE:['#f90','#fff'],LOW:['#48f','#fff'],MODIFIER:['#6b7494','#fff']};
function badge(imp){ const [bg,fg]=IC[imp]||['#ccc','#333']; return `<span class="badge" style="background:${bg};color:${fg}">${imp||'—'}</span>`; }
function esc(s){ if(s==null) return ''; return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }

// ── position index ────────────────────────────────────────────────────────
function buildPosIndex(data){
  const F={};  data.fields.forEach((f,i)=>F[f]=i);
  const idx={};
  data.rows.forEach((r,i)=>{
    const c=r[F.chrom]||'';  const p=parseInt(r[F.pos])||0;
    (idx[c]=idx[c]||[]).push([p,i]);
  });
  Object.values(idx).forEach(a=>a.sort((a,b)=>a[0]-b[0]));
  return idx;
}

function bisectLeft(arr, target){
  let lo=0, hi=arr.length;
  while(lo<hi){ const mid=(lo+hi)>>1; if(arr[mid][0]<target) lo=mid+1; else hi=mid; }
  return lo;
}

function domChrom(data){
  const F={};  data.fields.forEach((f,i)=>F[f]=i);
  const cc={};
  data.rows.forEach(r=>{ const c=r[F.chrom]||''; cc[c]=(cc[c]||0)+1; });
  return Object.entries(cc).sort((a,b)=>b[1]-a[1])[0]?.[0]||'';
}

// ── CART & EXPORT ─────────────────────────────────────────────────────────
function tagVariant(jsid, rowIdx, state) {
    const rowId = 'vrow-' + jsid + '-' + rowIdx;
    const row = document.getElementById(rowId);
    let cart = JSON.parse(localStorage.getItem('cart_' + jsid) || '{}');
    if (cart[rowId] === state) {
        delete cart[rowId];
        if(row) row.classList.remove('saved', 'discarded');
    } else {
        cart[rowId] = state;
        if(row) {
            row.classList.remove('saved', 'discarded');
            row.classList.add(state);
        }
    }
    localStorage.setItem('cart_' + jsid, JSON.stringify(cart));
}

function exportCartCSV(jsid) {
    const data = _varData[jsid];
    if (!data) return alert("Please load variants first.");
    const cart = JSON.parse(localStorage.getItem('cart_' + jsid) || '{}');
    let csv = "Triage_State," + data.fields.join(",") + "\\n";
    let count = 0;
    data.rows.forEach((r, i) => {
        const state = cart['vrow-' + jsid + '-' + i];
        if (state && state !== 'discarded' && state !== 'unreviewed') {
            const cleanRow = r.map(v => `"${(v||'').toString().replace(/"/g, '""')}"`).join(",");
            csv += `${state},${cleanRow}\\n`;
            count++;
        }
    });
    if(count === 0) return alert("No candidates saved in the cart yet!");
    const blob = new Blob([csv], {type: 'text/csv'});
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url; a.download = `triage_candidates_${jsid}.csv`;
    a.click(); URL.revokeObjectURL(url);
}

// ── row building ──────────────────────────────────────────────────────────
function makeRow(data, rowIdx, gi, F, jsid){
  const r = data.rows[rowIdx];
  const imp = r[F.impact] || '';
  const lof = r[F.lof] ? '<span class="lof-tag">LOF</span>' : '';
  const nmd = r[F.nmd] ? '<span class="nmd-tag">NMD</span>' : '';
  
  // Format sample columns to explicitly show reads
  const gtC = gi.map(ci => {
      const val = String(r[ci] || './. (?x)');
      const match = val.match(/^(.+?)\s*\((.*?)\)$/);
      if (match) {
          return `<td class="gt" style="white-space:nowrap;"><strong>${esc(match[1])}</strong> <span style="color:var(--muted); font-size:10px; margin-left:6px; background:rgba(255,255,255,0.05); padding:2px 4px; border-radius:3px;">Reads: ${esc(match[2].replace('x',''))}</span></td>`;
      }
      return `<td class="gt">${esc(val)}</td>`;
  }).join('');
  
  const tr = document.createElement('tr');
  const rowId = `vrow-${jsid}-${rowIdx}`;
  
  const cart = JSON.parse(localStorage.getItem('cart_'+jsid) || '{}');
  const initialClass = cart[rowId] ? ' ' + cart[rowId] : '';
  
  tr.className = 'var-row' + initialClass;
  tr.id = rowId;
  tr.dataset.impact = imp;
  tr.dataset.pos = r[F.pos] || '';
  tr.dataset.chrom = r[F.chrom] || '';
  tr.dataset.text = [r[F.vep_symbol], r[F.vep_consequence], r[F.vep_hgvsc],
                     r[F.ann_gene], r[F.ann_annotation], r[F.chrom]].join(' ').toLowerCase();

  let mutHits = 0, mutTotal = 0, sibHits = 0, sibTotal = 0;
  data.fields.forEach((colName, idx) => {
      const lowerCol = colName.toLowerCase();
      if (lowerCol.startsWith('s:')) {
          if (lowerCol.includes('mut')) { mutTotal++; if (r[idx] && r[idx].includes('1/1')) mutHits++; }
          if (lowerCol.includes('sib') || lowerCol.includes('wt')) { sibTotal++; if (r[idx] && r[idx].includes('1/1')) sibHits++; }
      }
  });
  let segScore = mutTotal > 0 ? (mutHits / mutTotal) - (sibTotal > 0 ? (sibHits / sibTotal) : 0) : 0;
  r._segScore = segScore;
  
  let segBadge = `<span class="badge" style="background:#6b7494;color:#fff">${segScore.toFixed(2)}</span>`;
  if (segScore >= 0.8) segBadge = `<span class="badge" style="background:#2ecc71;color:#fff">${segScore.toFixed(2)}</span>`;
  else if (segScore > 0.4) segBadge = `<span class="badge" style="background:#f39c12;color:#fff">${segScore.toFixed(2)}</span>`;

  const triageBtns = `
    <td style="text-align:center; white-space:nowrap; border-right: 1px solid var(--brd);" onclick="event.stopPropagation();">
      <button onclick="tagVariant('${jsid}', ${rowIdx}, 'saved')" style="background:#2ecc71; color:white; border:none; border-radius:3px; padding:2px 6px; cursor:pointer; font-weight:bold; margin-right:2px;" title="Save">✓</button>
      <button onclick="tagVariant('${jsid}', ${rowIdx}, 'discarded')" style="background:#e74c3c; color:white; border:none; border-radius:3px; padding:2px 6px; cursor:pointer; font-weight:bold;" title="Discard">✕</button>
    </td>
  `;

  tr.innerHTML =
    triageBtns +
    `<td style="text-align:center;">${segBadge}</td>` +
    `<td>${badge(imp)}${lof}${nmd}</td>` +
    `<td>${esc(r[F.chrom])}</td>` +
    `<td>${esc(r[F.pos])}</td>` +
    `<td>${esc(r[F.ref])}</td>` +
    `<td>${esc(r[F.alt])}</td>` +
    `<td><strong>${esc(r[F.info_dp] || '0')}</strong></td>` +
    `<td>${esc(r[F.qual])}</td>` +
    `<td>${esc(r[F.filter])}</td>` +
    `<td>${esc(r[F.vep_symbol] || '—')}</td>` +
    `<td class="sc">${esc(r[F.vep_consequence] || '—')}</td>` +
    `<td class="sc">${esc(r[F.vep_hgvsc] || '—')}</td>` +
    `<td class="sc">${esc(r[F.vep_hgvsp] || '—')}</td>` +
    `<td>${esc(r[F.vep_aa] || '—')}</td>` +
    `<td class="tc">${r[F.vep_n] || '—'}</td>` +
    `<td>${esc(r[F.ann_gene] || '—')}</td>` +
    `<td class="sc">${esc(r[F.ann_annotation] || '—')}</td>` +
    `<td class="tc">${r[F.ann_n] || '—'}</td>` +
    gtC;
    
  tr.addEventListener('click',()=>navigateToCoordinate(jsid, r[F.chrom], parseInt(r[F.pos]) || 0));
  tr.addEventListener('mouseover',()=>updateGffMarker(jsid, r[F.chrom], parseInt(r[F.pos]) || 0));
  return tr;
}

function buildHeader(data, headId, sampleNames, jsid){
  const F={}; data.fields.forEach((f,i)=>F[f]=i);
  
  const tooltips = {
      seg: "Continuous score calculating fidelity of mutants against sibling noise.",
      impact: "Predicted molecular variant severity impact classification.",
      chrom: "Chromosome or Contig identifier.",
      pos: "1-based genomic coordinate position.",
      dp: "Total sequence reads covering this variant locus site. Note: This number will not match the absolute read depth per sample. Broadly, it is a aggregate score based on multiple variables. However, it can be informative as one (of many) confidence scores.",
      qual: "Phred-scaled base call quality configuration score.",
      filter: "Filter status validation rule flags.",
      gene: "Ensembl Gene Symbol mapping.",
      cons: "Primary VEP structural consequence array map.",
      ann_gene: "SnpEff Gene Symbol.",
      ann_cons: "Primary SnpEff annotation."
  };

  const makeTh = (label, sortKey, type, tooltip) => {
      const sortAction = sortKey ? ` onclick="sortTable('${jsid}', '${sortKey}', '${type}')" style="cursor:pointer; text-decoration:underline; text-decoration-style:dotted; text-underline-offset:3px;"` : '';
      const ttHtml = tooltip ? ` <span class="tt-icon" data-tt="${tooltip}">?</span>` : '';
      return `<th><div class="th-wrap"><span${sortAction}>${label} ${sortKey ? '↕' : ''}</span>${ttHtml}</div></th>`;
  };

  const gtH = sampleNames.map((s, idx) => {
      const fieldIdx = F[`s:${s}`];
      return makeTh(`${esc(s)}`, `gt_${fieldIdx}`, 'gt', `Genotype call and local read depth for sample ${esc(s)}. This should correspond more closely to the true number of reads mapping to the region.`);
  }).join('');

  document.getElementById(headId).innerHTML = `<tr>` +
    `<th style="width:80px; text-align:center;"><div class="th-wrap" style="justify-content:center;"><span>Triage</span><span class="tt-icon" data-tt="Tag variants for curation. Preserved on exports.">?</span></div></th>` +
    makeTh('Seg', 'segScore', 'num', tooltips.seg) +
    makeTh('Impact', 'impact', 'impact', tooltips.impact) +
    makeTh('Chrom', 'chrom', 'str', tooltips.chrom) +
    makeTh('Pos', 'pos', 'num', tooltips.pos) +
    `<th>Ref</th><th>Alt</th>` +
    makeTh('Depth', 'info_dp', 'num', tooltips.dp) +
    makeTh('Qual', 'qual', 'num', tooltips.qual) +
    `<th>Filter</th>` +
    makeTh('VEP Sym', 'vep_symbol', 'str', tooltips.gene) +
    makeTh('Consequence', 'vep_consequence', 'str', tooltips.cons) +
    `<th>HGVSc</th><th>HGVSp</th><th>AA</th><th>VEP tx</th>` +
    makeTh('ANN Gene', 'ann_gene', 'str', tooltips.ann_gene) +
    makeTh('Annotation', 'ann_annotation', 'str', tooltips.ann_cons) +
    `<th>ANN tx</th>` + gtH + `</tr>`;
}

// ── Core Array Sorting Engine ──────────────────────────────────────────────
function sortTable(jsid, fieldKey, type) {
    const data = _varData[jsid];
    let activeIndices = _currentIndices[jsid];
    if (!data || !activeIndices) return;

    _sortState[jsid] = _sortState[jsid] || { key: '', order: 'asc' };
    const state = _sortState[jsid];
    state.order = (state.key === fieldKey && state.order === 'asc') ? 'desc' : 'asc';
    state.key = fieldKey;
    const dir = state.order === 'asc' ? 1 : -1;

    const F = {}; data.fields.forEach((f, i) => F[f] = i);
    const impactVal = {"HIGH":4, "MODERATE":3, "LOW":2, "MODIFIER":1, "":0};

    activeIndices.sort((idxA, idxB) => {
        const a = data.rows[idxA];
        const b = data.rows[idxB];
        let valA, valB;

        if (fieldKey === 'segScore') {
            valA = a._segScore || 0; valB = b._segScore || 0;
        } else if (type === 'gt') {
            const colIdx = parseInt(fieldKey.split('_')[1]);
            valA = parseInt((String(a[colIdx]||'').match(/\((\d+)x\)/) || [0,0])[1]);
            valB = parseInt((String(b[colIdx]||'').match(/\((\d+)x\)/) || [0,0])[1]);
        } else {
            valA = a[F[fieldKey]]; valB = b[F[fieldKey]];
        }

        if (type === 'num' || type === 'gt') {
            return ((parseFloat(valA) || 0) - (parseFloat(valB) || 0)) * dir;
        } else if (type === 'impact') {
            return ((impactVal[valA||""]||0) - (impactVal[valB||""]||0)) * dir;
        } else {
            const strA = String(valA||"").toLowerCase();
            const strB = String(valB||"").toLowerCase();
            return strA < strB ? -1 * dir : (strA > strB ? 1 * dir : 0);
        }
    });

    const gi = (window._snames[jsid] || []).map(s => F[`s:${s}`]);
    const tbody = document.getElementById(`tb-${jsid}`);
    if (tbody) {
        tbody.innerHTML = '';
        const frag = document.createDocumentFragment();
        activeIndices.forEach(ri => frag.appendChild(makeRow(data, ri, gi, F, jsid)));
        tbody.appendChild(frag);
        runStratification(jsid); // Re-apply opacities to sorted rows
    }
}

// ── Stratification Visualizer ──────────────────────────────────────────────
function runStratification(jsid) {
    const modeEl = document.getElementById('preset-mode-' + jsid);
    const mode = modeEl ? modeEl.value : 'relaxed';
    const minDp = parseInt((document.getElementById('min-dp-' + jsid) || {}).value || '0', 10);
    const data = _varData[jsid]; if(!data) return;
    const F = {}; data.fields.forEach((f,i) => F[f] = i);

    const mutCols = [];
    data.fields.forEach((colName, idx) => { if (colName.toLowerCase().includes('mut')) mutCols.push(idx); });

    document.querySelectorAll(`#tw-${jsid} .var-row`).forEach(tr => {
        // Skip explicitly triaged rows, force them to stay fully opaque
        if (tr.classList.contains('saved') || tr.classList.contains('discarded')) return;
        
        const rIdx = parseInt(tr.id.split('-').pop(), 10);
        const r = data.rows[rIdx]; if(!r) return;
        
        let score = 0;
        const imp = r[F.impact] || '';
        if (imp === 'HIGH') score += 2;
        if (imp === 'MODERATE') score += 1;
        if (r._segScore >= 0.8) score += 2;

        let passDepth = true;
        if (mutCols.length > 0 && minDp > 0) {
            mutCols.forEach(idx => {
                const match = String(r[idx]||'').match(/\((\d+)x\)/);
                if (match && parseInt(match[1], 10) < minDp) passDepth = false;
            });
        }

        // Clean previous tiers
        tr.classList.remove('tier-1', 'tier-2', 'tier-3');
        
        // Force assignment
        if (mode === 'strict') {
            if (score >= 3 && passDepth) tr.classList.add('tier-1');
            else tr.classList.add('tier-3'); 
        } 
        else if (mode === 'segregation') {
            if (r._segScore >= 0.8 && passDepth) tr.classList.add('tier-1');
            else if (r._segScore > 0) tr.classList.add('tier-2');
            else tr.classList.add('tier-3');
        } 
        else {
            if (!passDepth) tr.classList.add('tier-3');
            else if (score >= 3) tr.classList.add('tier-1');
            else if (score > 0) tr.classList.add('tier-2');
            else tr.classList.add('tier-3');
        }
    });
}


// ── windowed render (±1000 rows around a position, by sorted index) ─────────
function renderWindow(jsid, chrom, centerPos, sampleNames){
  const data=_varData[jsid];  if(!data) return;
  const idx=_posIndex[jsid];  if(!idx) return;
  const F={};  data.fields.forEach((f,i)=>F[f]=i);
  const gi=sampleNames.map(s=>F[`s:${s}`]);
  const arr=idx[chrom]||[];
  const ci=bisectLeft(arr,centerPos);
  const lo=Math.max(0,ci-1000);
  const hi=Math.min(arr.length,ci+1000);
  const rowIndices=arr.slice(lo,hi).map(([,ri])=>ri);
  _viewMode[jsid]='window';
  setStatus(jsid,`Showing ${rowIndices.length} variants around ${chrom}:${centerPos.toLocaleString()} (±1000)`);
  const tbody=document.getElementById('tb-'+jsid);
  tbody.innerHTML='';
  const frag=document.createDocumentFragment();
  rowIndices.forEach(ri=>frag.appendChild(makeRow(data,ri,gi,F,jsid)));
  tbody.appendChild(frag);
  document.getElementById('tw-'+jsid).style.display='';
  _currentIndices[jsid] = rowIndices;
  runStratification(jsid);
}

// ── full chunked render (SNAP / show all) ──────────────────────────────────
function renderAll(jsid, filterImpactVal, sampleNames){
  const data=_varData[jsid];  if(!data) return;
  const F={};  data.fields.forEach((f,i)=>F[f]=i);
  const gi=sampleNames.map(s=>F[`s:${s}`]);
  _viewMode[jsid]='all';

  let indices=[];
  data.rows.forEach((r,i)=>{
    if(!filterImpactVal||filterImpactVal==='ALL'||r[F.impact]===filterImpactVal) indices.push(i);
  });

  const tbody=document.getElementById('tb-'+jsid);
  tbody.innerHTML='';
  document.getElementById('vstatus-'+jsid).style.display='';
  document.getElementById('tw-'+jsid).style.display='';

  const total=indices.length;
  let i=0;
  function nextChunk(){
    const chunkSize=Math.min(200+Math.floor(i/500)*100, 2000);
    const frag=document.createDocumentFragment();
    const end=Math.min(i+chunkSize,total);
    for(;i<end;i++) frag.appendChild(makeRow(data,indices[i],gi,F,jsid));
    tbody.appendChild(frag);
    const pct=total?Math.round((i/total)*100):100;
    const pbar=document.getElementById('pbar-'+jsid);
    const ppct=document.getElementById('ppct-'+jsid);
    if(pbar) pbar.style.width=pct+'%';
    if(ppct) ppct.textContent=pct+'%';
    setStatus(jsid,`Rendering… ${i.toLocaleString()} / ${total.toLocaleString()}`);
    if(i<total){ setTimeout(nextChunk,0); }
    else{ document.getElementById('vstatus-'+jsid).style.display='none'; _currentIndices[jsid] = indices; applyFilters('tw-'+jsid,jsid); runStratification(jsid); }
  }
  nextChunk();
}

function setStatus(jsid, msg){
  const el=document.getElementById('vmsg-'+jsid);
  if(el) el.textContent=msg;
}

// ── filter / search (works in both modes) ─────────────────────────────────
const activeFilters={};
function filterImpact(btn, impact, tblId, jsid){
  activeFilters[tblId]=impact;
  const fbar=document.getElementById('fbar-'+jsid);
  if(fbar) fbar.querySelectorAll('.f-btn').forEach(b=>b.classList.remove('active'));
  btn.classList.add('active');
  applyFilters(tblId,jsid);
}
function applyFilters(tblId, jsid){
  const wrap=document.getElementById(tblId);
  if(!wrap) return;
  const impact=activeFilters[tblId]||'ALL';
  const search=(document.getElementById('sterm-'+jsid)||{value:''}).value.toLowerCase();
  const posRaw=(document.getElementById('spos-'+jsid)||{value:''}).value.trim();
  let posMin=null, posMax=null;
  if(posRaw){
    const clean=posRaw.replace(/^.*?:/,'');
    if(clean.includes('-')){ const p=clean.split('-'); posMin=parseInt(p[0]); posMax=parseInt(p[1]); }
    else if(clean.includes('..')){ const p=clean.split('..'); posMin=parseInt(p[0]); posMax=parseInt(p[1]); }
    else{ posMin=posMax=parseInt(clean); }
  }
  wrap.querySelectorAll('.var-row').forEach(row=>{
    const impOk=impact==='ALL'||row.dataset.impact===impact;
    const searchOk=!search||(row.dataset.text||'').includes(search);
    let posOk=true;
    if(posMin!==null&&!isNaN(posMin)){
      const rp=parseInt(row.dataset.pos);
      posOk=!isNaN(rp)&&rp>=posMin&&(posMax===null||isNaN(posMax)||rp<=posMax);
    }
    row.style.display=(impOk&&searchOk&&posOk)?'':'none';
  });
}

// ── IGV Integration ───────────────────────────────────────────────────────
function renderIGV(jsid, configStr) {
    if(!configStr) return;
    const igvDiv = document.getElementById('igv-container-' + jsid);
    if(!igvDiv || typeof igv === 'undefined') return;
    const options = JSON.parse(configStr);
    igv.createBrowser(igvDiv, options);
}

// ── GFF track ─────────────────────────────────────────────────────────────
const FEAT_COLORS={gene:'#4f8ef7',mRNA:'#7c5be8',CDS:'#2a9d8f',exon:'#e76f51',transcript:'#7c5be8',pseudogene:'#888',ncRNA:'#f4a261'};
const IMPACT_COLORS={HIGH:'#f44',MODERATE:'#f90',LOW:'#48f',MODIFIER:'#6b7494'};
const IMPACT_STYLE={HIGH:{y1:25,w:3.5},MODERATE:{y1:40,w:2.5},LOW:{y1:50,w:2},MODIFIER:{y1:55,w:1.5}};

function renderGffTrack(features, data, jsid){
  const svg=document.getElementById('gff-svg-'+jsid);
  const coords=document.getElementById('gff-coords-'+jsid);
  if(!svg) return;
  const F={};  data.fields.forEach((f,i)=>F[f]=i);
  const positions=[];
  for(let i=0;i<data.rows.length;i++){ const p=parseInt(data.rows[i][F.pos])||0; if(p>0) positions.push(p); }
  if(!positions.length) return;
  let minPos=positions[0], maxPos=positions[0];
  for(let i=1;i<positions.length;i++){ if(positions[i]<minPos) minPos=positions[i]; if(positions[i]>maxPos) maxPos=positions[i]; }
  const span=maxPos-minPos||1;
  const dc=domChrom(data);
  const norm=c=>String(c).toLowerCase().replace(/^chr/,'');
  const dcn=norm(dc);
  const relevant=features.filter(f=>norm(f.chrom)===dcn);
  coords.textContent=`${dc}:${minPos.toLocaleString()}–${maxPos.toLocaleString()}`;
  const W=svg.getBoundingClientRect().width||800;
  const toX=pos=>((pos-minPos)/span)*(W-20)+10;
  const toPos=x=>Math.round(((x-10)/(W-20))*span+minPos);
  svg.innerHTML='';
  const ns='http://www.w3.org/2000/svg';
  // axis
  const axis=document.createElementNS(ns,'line');
  axis.setAttribute('x1',10);axis.setAttribute('y1',70);axis.setAttribute('x2',W-10);axis.setAttribute('y2',70);
  axis.setAttribute('stroke','#2a2f45');axis.setAttribute('stroke-width','1');
  svg.appendChild(axis);
  // variant lines
  data.rows.forEach((r,i)=>{
    const imp=r[F.impact]||'MODIFIER';
    const pos=parseInt(r[F.pos])||0;  if(!pos) return;
    const st=IMPACT_STYLE[imp]||IMPACT_STYLE.MODIFIER;
    const line=document.createElementNS(ns,'line');
    line.setAttribute('x1',toX(pos));line.setAttribute('x2',toX(pos));
    line.setAttribute('y1',st.y1);line.setAttribute('y2',70);
    line.setAttribute('stroke',IMPACT_COLORS[imp]||'#6b7494');
    line.setAttribute('stroke-width',st.w);
    line.style.cursor='pointer';
    line.addEventListener('click',()=>navigateToCoordinate(jsid,r[F.chrom],pos,null));
    svg.appendChild(line);
  });
  // GFF features
  const ROW_Y={gene:15,mRNA:32,CDS:49,exon:49,transcript:32};
  const ROW_H={gene:10,mRNA:8,CDS:10,exon:6,transcript:8};
  relevant.forEach(feat=>{
    const x1=toX(feat.start),x2=toX(feat.end),w=Math.max(2,x2-x1);
    const ry=ROW_Y[feat.type]??49,rh=ROW_H[feat.type]??8;
    const rect=document.createElementNS(ns,'rect');
    rect.setAttribute('x',x1);rect.setAttribute('y',ry);rect.setAttribute('width',w);rect.setAttribute('height',rh);
    rect.setAttribute('fill',FEAT_COLORS[feat.type]||'#6b7494');rect.setAttribute('opacity','0.85');rect.setAttribute('rx','2');
    rect.style.cursor='pointer';
    const title=document.createElementNS(ns,'title');
    title.textContent=`${feat.name||'?'} [${feat.type}] ${feat.chrom}:${feat.start}-${feat.end} (${feat.strand}) — click to navigate`;
    rect.appendChild(title);
    rect.addEventListener('click',()=>navigateToCoordinate(jsid,dc,feat.start,feat.end));
    svg.appendChild(rect);
    if(w>40&&feat.name&&feat.type==='gene'){
      const txt=document.createElementNS(ns,'text');
      txt.setAttribute('x',x1+2);txt.setAttribute('y',ry+rh-1);txt.setAttribute('font-size','9');txt.setAttribute('fill','#d0d6f0');txt.setAttribute('pointer-events','none');
      txt.textContent=feat.name.length>12?feat.name.slice(0,11)+'…':feat.name;
      svg.appendChild(txt);
    }
  });
  // click on empty axis space → navigate to that position
  svg.addEventListener('click',e=>{
    if(e.target!==svg&&e.target.tagName!=='line') return;
    const x=e.clientX-svg.getBoundingClientRect().left;
    const pos=toPos(x);
    if(pos>=minPos&&pos<=maxPos) navigateToCoordinate(jsid,dc,pos,null);
  });
  // hover marker
  const marker=document.createElementNS(ns,'line');
  marker.id='gff-marker-'+jsid;
  marker.setAttribute('x1',0);marker.setAttribute('y1',0);marker.setAttribute('x2',0);marker.setAttribute('y2',75);
  marker.setAttribute('stroke','#ff4444');marker.setAttribute('stroke-width','1.5');marker.setAttribute('stroke-dasharray','3,2');
  marker.style.display='none';
  svg.appendChild(marker);
  svg._gffScale={minPos,span,W,dc};
}

function updateGffMarker(jsid,chrom,pos){
  const svg=document.getElementById('gff-svg-'+jsid);
  if(!svg||!svg._gffScale) return;
  const {minPos,span,W}=svg._gffScale;
  const marker=document.getElementById('gff-marker-'+jsid);
  if(!marker) return;
  const x=((pos-minPos)/span)*(W-20)+10;
  marker.setAttribute('x1',x);marker.setAttribute('x2',x);marker.style.display='';
  const coords=document.getElementById('gff-coords-'+jsid);
  if(coords) coords.textContent=`${chrom}:${pos.toLocaleString()} (hover)`;
}

// Navigate to coordinate: windowed render + marker + position filter + igv sync
function navigateToCoordinate(jsid, chrom, startPos, endPos, sampleNames){
  const posInp=document.getElementById('spos-'+jsid);
  if(posInp) posInp.value=endPos?`${chrom}:${startPos}-${endPos}`:`${chrom}:${startPos}`;
  
  const sn=sampleNames||window._snames[jsid]||[];
  renderWindow(jsid,chrom,startPos,sn);
  updateGffMarker(jsid,chrom,startPos);
  
  // Sync back to IGV Frame (Pad by 150bp for visual context)
  if (window.igvBrowsers && window.igvBrowsers[jsid]) {
      window._isSyncing = true;
      const targetStart = Math.max(1, startPos - 150);
      const targetEnd = startPos + 150;
      window.igvBrowsers[jsid].search(`${chrom}:${targetStart}-${targetEnd}`)
          .finally(() => { window._isSyncing = false; });
  }

  setTimeout(()=>{
    const first=document.querySelector('#tw-'+jsid+' .var-row:not([style*="none"])');
    if(first) first.scrollIntoView({behavior:'smooth',block:'center'});
  },50);
}

// Web Worker parse helper — inline worker via Blob so no separate .js file needed
function parseJsonInWorker(text){
  return new Promise((resolve,reject)=>{
    const src=`self.onmessage=e=>{ try{ postMessage(JSON.parse(e.data)); }catch(err){ postMessage({__error:err.message}); } };`;
    const blob=new Blob([src],{type:'application/javascript'});
    const url=URL.createObjectURL(blob);
    const worker=new Worker(url);
    worker.onmessage=e=>{ URL.revokeObjectURL(url); worker.terminate(); if(e.data&&e.data.__error) reject(new Error(e.data.__error)); else resolve(e.data); };
    worker.onerror=e=>{ URL.revokeObjectURL(url); worker.terminate(); reject(e); };
    worker.postMessage(text);
  });
}
"""

# ── HTML rendering ────────────────────────────────────────────────────────────

def _gallery_html(images, title, section_id):
    body = ""
    if not images:
        body = '<p class="muted">No images found.</p>'
    else:
        figs = []
        for img in images:
            try:
                src  = encode_image(img)
                figs.append(f'<figure class="img-figure" onclick="openLightbox(this)">'
                             f'<img src="{src}" alt="{_escape(img.name)}" loading="lazy">'
                             f'<figcaption>{_escape(img.name)}</figcaption></figure>')
            except Exception as e:
                figs.append(f'<p class="muted">[Could not load {img.name}: {e}]</p>')
        body = f'<div class="gallery">{"".join(figs)}</div>'
    count = f'<span class="count">({len(images)})</span>' if images else ''
    return (f'<div class="section">'
            f'<div class="section-hdr" onclick="toggleSection(\'{section_id}\')">'
            f'<h3>{title} {count}</h3>'
            f'<span class="section-arrow" id="ar-{section_id}">▶</span></div>'
            f'<div class="section-body" id="sb-{section_id}">{body}</div></div>')

# ── GFF extraction ────────────────────────────────────────────────────────────

def extract_gff_features(gff_path: Path, chrom: str, start: int, end: int) -> list:
    if not gff_path or not gff_path.exists():
        return []
    import subprocess
    try:
        result = subprocess.run(
            ["tabix", str(gff_path), f"{chrom}:{start}-{end}"],
            capture_output=True, text=True, timeout=30
        )
        if result.returncode != 0:
            alt = chrom.lstrip("chr") if chrom.startswith("chr") else f"chr{chrom}"
            result = subprocess.run(
                ["tabix", str(gff_path), f"{alt}:{start}-{end}"],
                capture_output=True, text=True, timeout=30
            )
        features = []
        for line in result.stdout.splitlines():
            if not line or line.startswith("#"): continue
            parts = line.split("\t")
            if len(parts) < 9: continue
            fchrom, _src, ftype, fstart, fend, _score, strand = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5], parts[6]
            attrs = parts[8]
            name = ""
            for pat in [r'Name=([^;]+)', r'gene_name=([^;]+)', r'gene_id=([^;]+)', r'ID=([^;]+)']:
                m = re.search(pat, attrs)
                if m: name = m.group(1); break
            try:
                features.append({"chrom": fchrom, "start": int(fstart), "end": int(fend),
                                  "strand": strand, "type": ftype, "name": name})
            except ValueError:
                continue
        return features
    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        print(f"  [WARN] tabix GFF extraction failed: {e}", file=sys.stderr)
        return []

def get_vcf_coord_range(json_path: Path) -> dict:
    if not json_path.exists(): return {}
    try:
        with open(json_path, encoding="utf-8") as f:
            data = json.load(f)
        F = {name: i for i, name in enumerate(data["fields"])}
        fi_chrom, fi_pos = F.get("chrom"), F.get("pos")
        if fi_chrom is None or fi_pos is None: return {}
        ranges = {}
        for row in data["rows"]:
            chrom = row[fi_chrom]
            try: pos = int(row[fi_pos])
            except (ValueError, TypeError): continue
            if chrom not in ranges:
                ranges[chrom] = [pos, pos]
            else:
                if pos < ranges[chrom][0]: ranges[chrom][0] = pos
                if pos > ranges[chrom][1]: ranges[chrom][1] = pos
        return {c: (v[0], v[1]) for c, v in ranges.items()}
    except Exception as e:
        print(f"  [WARN] coord range read failed: {e}", file=sys.stderr)
        return {}

# ── IGV configuration ─────────────────────────────────────────────────────────

def get_bam_locus(bam_path: Path) -> str:
    import subprocess
    try:
        idx_res = subprocess.run(["samtools", "idxstats", str(bam_path)], capture_output=True, text=True)
        best_chrom, max_reads = "", 0
        for line in idx_res.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) >= 3 and int(parts[2]) > max_reads:
                max_reads, best_chrom = int(parts[2]), parts[0]
        if not best_chrom: return "all"
        view_res = subprocess.run(["samtools", "view", str(bam_path), best_chrom], capture_output=True, text=True)
        first_pos = int(view_res.stdout.splitlines()[0].split("\t")[3])
        return f"{best_chrom}:{max(1, first_pos - 2000)}-{first_pos + 2000}"
    except Exception: return "all"

def build_igv_config(sample_dir: Path, fasta_name: str = None, gff_name: str = None) -> str:
    tracks = []
    prefix = sample_dir.name
    
    go_dir = sample_dir / "genome_overview"
    if go_dir.is_dir():
        for f in sorted(go_dir.iterdir()):
            if f.suffix in [".tdf", ".bedgraph", ".bw"]:
                fmt = "tdf" if f.suffix == ".tdf" else "bedgraph"
                tracks.append({"name": f"Overview: {f.name}", "url": f"{prefix}/genome_overview/{f.name}", "type": "wig", "format": fmt, "order": 10})

    nt_dir = sample_dir / "nucleotide_view"
    if nt_dir.is_dir():
        for sub in ["mut", "wt", "sib"]:
            sub_dir = nt_dir / sub
            if sub_dir.is_dir():
                for f in sorted(sub_dir.iterdir()):
                    if f.suffix == ".bam":
                        order = 50 if sub == "mut" else 60
                        tracks.append({"name": f"Reads ({sub.upper()}): {f.name}", "url": f"{prefix}/nucleotide_view/{sub}/{f.name}", "indexURL": f"{prefix}/nucleotide_view/{sub}/{f.name}.bai", "format": "bam", "type": "alignment", "order": order})

    if gff_name:
        is_gz = gff_name.endswith('.gz')
        tracks.append({
            "name": "Custom Annotations",
            "type": "annotation",
            "format": "gtf" if "gtf" in gff_name.lower() else "gff3",
            "url": gff_name,
            "indexURL": f"{gff_name}.tbi" if is_gz else None,
            "displayMode": "EXPANDED",
            "color": "#3498db",
            "height": 120,
            "visibilityWindow": 1000000
        })

    if not tracks: return ""
    locus = "all"
    first_bam = next((sample_dir.parent / t["url"] for t in tracks if t["format"] == "bam"), None)
    if first_bam and first_bam.exists(): locus = get_bam_locus(first_bam)     
    
    config = {"locus": locus, "tracks": tracks}
    if fasta_name: config["reference"] = {"id": "danRer7_local", "name": "Zebrafish (danRer7)", "fastaURL": fasta_name, "indexURL": f"{fasta_name}.fai"}
    else: config["reference"] = {"id": "danRer7_ucsc", "name": "Zebrafish (danRer7)", "fastaURL": "https://hgdownload.soe.ucsc.edu/goldenPath/danRer7/bigZips/danRer7.fa.gz", "indexURL": "https://hgdownload.soe.ucsc.edu/goldenPath/danRer7/bigZips/danRer7.fa.gz.fai"}
    return json.dumps(config)

# ── VCF section HTML (inline on sample page, non-collapsible) ─────────────────

def _vcf_section_html(vcf_path, jsid, max_vcf_rows, fasta_name, gff_name, has_gff: bool):
    """Return HTML for one VCF block (non-collapsible, inline variant table)."""
    _, col_headers, sample_names, csq_fmt, ann_fmt = read_vcf_header(vcf_path)
    no_vep = '' if csq_fmt else '<p class="muted" style="margin:.2rem 0;font-size:11px">No VEP/CSQ annotations in this VCF.</p>'
    no_ann = '' if ann_fmt else '<p class="muted" style="margin:.2rem 0;font-size:11px">No SnpEff/ANN annotations in this VCF.</p>'
    if not col_headers:
        return f'<div class="vcf-block"><h4 class="vcf-name">{_escape(vcf_path.name)}</h4><p class="muted">Could not read VCF header.</p></div>'

    gt_header_js = json.dumps(sample_names)
    trunc_note   = (f'<p class="muted" style="display:none" id="tn-{jsid}">'
                    f'Showing first {max_vcf_rows} variants.</p>') if max_vcf_rows else ''
    tbl_wrap_id  = f"tw-{jsid}"

    sample_dir = vcf_path.parent
    igv_config_json = build_igv_config(sample_dir, fasta_name, gff_name)
    igv_html = f"""
    <div id="igv-container-{jsid}" style="margin: 1rem 0; border: 1px solid var(--brd); border-radius: var(--r); background: #ffffff; color: #111111; padding: 5px;"></div>
    <script>
      document.addEventListener("DOMContentLoaded", function() {{
          if (typeof igv !== 'undefined') {{
              window.igvBrowsers = window.igvBrowsers || {{}};
              igv.createBrowser(document.getElementById('igv-container-{jsid}'), {igv_config_json})
                 .then(function(b) {{
                     window.igvBrowsers['{jsid}'] = b;
                     // Sync IGV locus panning back to the custom GFF readout
                     b.on('locuschange', function(ref) {{
                         if (window._isSyncing) return;
                         window._isSyncing = true;
                         const f = Array.isArray(ref) ? ref[0] : ref;
                         if (f) {{
                             const c = document.getElementById('gff-coords-{jsid}');
                             if (c) c.textContent = `${{f.chr}}:${{Math.floor(f.start)}}-${{Math.ceil(f.end)}}`;
                         }}
                         window._isSyncing = false;
                     }});
                 }});
          }}
      }});
    </script>
    """ if igv_config_json else ""

    gff_track_html = f"""
    <div id="gff-wrap-{jsid}" style="display:none; margin: .5rem 0;">
      <div class="gff-sticky" id="gff-sticky-{jsid}" style="border: 1px solid var(--brd); border-radius: var(--r); background: var(--surf); padding: .75rem;">
        <div class="gff-header" style="display: flex; justify-content: space-between; font-size: 11px; margin-bottom: .4rem;">
          <div class="th-wrap">
            <span class="gff-title" style="font-weight: bold; color: var(--acc);">INTERACTIVE GFF TRACK</span>
            <span class="tt-icon" data-tt="Interactive coordinate map. Click coloured blocks or ticks to jump the view to that locus and populate the table below with variants in that region. Hover over upper blue boxes for gene details. Red/orange/lower blue/lower grey regions correspond to predicted variant impacts of High/Moderate/Low/Modifier respectively. Note that two separate variant prediction algorithms are used, and while variant impact predictions are often in concordance, we display the highest impact call in cases of disagreement.">?</span>
          </div>
          <span class="gff-coords" id="gff-coords-{jsid}" style="font-family: var(--mono); color: var(--txt-muted);"></span>
        </div>
        <svg id="gff-svg-{jsid}" class="gff-svg" width="100%" height="80" style="background: var(--bg); display: block; overflow: visible;"></svg>
      </div>
    </div>""" if has_gff else ''

    return f"""
<div class="vcf-block">
  <h4 class="vcf-name">{_escape(vcf_path.name)}</h4>
  {no_vep}{no_ann}
  {igv_html}
  {gff_track_html}
  {trunc_note}
  
  <div style="margin-bottom:.75rem; display:flex; gap:.5rem;">
    <button class="load-btn" onclick="exportCartCSV('{jsid}')" style="background: rgba(46, 204, 113, 0.15); border-color: #2ecc71; color: #2ecc71;">
      Export Curated CSV
    </button>
  </div>

<div class="filter-bar" id="fbar-{jsid}" style="display: flex; flex-direction: column; gap: 10px;">
    <div style="display:flex; align-items:center; gap: 15px; flex-wrap:wrap;">
        <div style="display: inline-flex; align-items: center; gap: 6px; font-weight: bold; font-size: 12px; color: var(--txt);">
            <span>Preset Mode</span>
            <span class="tt-icon" data-tt="Toggles strictness tiers. 'Strict' requires high impact & passing depth. 'Segregation' prioritizes mutant vs. sibling ratio.">?</span>
            <select id="preset-mode-{jsid}" class="search-box" onchange="runStratification('{jsid}')" style="margin-left: 4px; padding: 2px 4px; width: auto;">
                <option value="relaxed">Relaxed</option>
                <option value="strict">Strict</option>
                <option value="segregation">Segregation</option>
            </select>
        </div>

        <div style="display: inline-flex; align-items: center; gap: 6px; font-weight: bold; font-size: 12px; color: var(--txt);">
            <span>Min Depth</span>
            <span class="tt-icon" data-tt="Hides variants where the Local Read Depth for the mutant samples falls below this threshold.">?</span>
            <input type="number" id="min-dp-{jsid}" value="0" min="0" class="search-box" onchange="runStratification('{jsid}')" style="margin-left: 4px; width: 60px; padding: 2px;">
        </div>
    </div>

    <div style="display:flex; align-items:center; gap: 10px; flex-wrap:wrap;">
        <span style="display: inline-flex; align-items: center; gap: 6px; font-weight: bold; font-size: 12px; color: var(--txt);">
            Filter
            <span class="tt-icon" data-tt="Regex/Text search. Filters by Gene Name, Consequence, or Chromosome.">?</span>
        </span>
        <button class="f-btn active" onclick="filterImpact(this,'ALL','{tbl_wrap_id}','{jsid}')">All</button>
        <button class="f-btn" onclick="filterImpact(this,'HIGH','{tbl_wrap_id}','{jsid}')">HIGH</button>
        <button class="f-btn" onclick="filterImpact(this,'MODERATE','{tbl_wrap_id}','{jsid}')">MODERATE</button>
        <button class="f-btn" onclick="filterImpact(this,'LOW','{tbl_wrap_id}','{jsid}')">LOW</button>
        <button class="f-btn" onclick="filterImpact(this,'MODIFIER','{tbl_wrap_id}','{jsid}')">MODIFIER</button>
        <input class="search-box" id="sterm-{jsid}" type="text" placeholder="Search gene / consequence…" oninput="applyFilters('{tbl_wrap_id}','{jsid}'); window.runStratification && runStratification('{jsid}');" style="width:200px">
        <input class="search-box" id="spos-{jsid}" type="text" placeholder="Pos range e.g. chr13:335-337" oninput="applyFilters('{tbl_wrap_id}','{jsid}')" style="width:200px">
    </div>
  </div>
  <p>
    The coordinate and search gene functionality work <b><i>within the scope of the selected interactive region only</i></b>. Click column headers to sort, and mouse over question marks for helptext. Opacity is automatically imposed on low-scoring variants, but all are shown for transparency. You can "freeze" specific variants of interest in the triage tab and export a file compatible with excel, if desired. Selecting a region on either the IGV browser or Interactive GFF Track will update the table, and propagate the coordinate of interest across other tracks. Click <b>Load variants</b> below to proceed. (Variants are not automatically loaded as this can be a heavy operation.)
  </p>
  
  <div id="vload-{jsid}" style="margin-bottom:.75rem">
    <button class="load-btn" onclick="loadVcf_{jsid}()">▶ Load variants</button>
    <p class="muted" style="margin-top:.4rem;font-size:11px">
    </p>
  </div>
  <div id="vstatus-{jsid}" style="display:none;padding:.5rem 0">
    <span class="muted" id="vmsg-{jsid}">Loading…</span>
  </div>
  <div id="{tbl_wrap_id}" style="display:none">
    <div class="tbl-scroll">
    <table class="var-table">
      <thead id="th-{jsid}"></thead>
      <tbody id="tb-{jsid}"></tbody>
    </table>
    </div>
  </div>
  <script>
  (function(){{
    const SNAMES  = {gt_header_js};
    const HAS_GFF = {'true' if has_gff else 'false'};
    const JSID    = '{jsid}';
    window._snames = window._snames || {{}};
    window._snames[JSID] = SNAMES;

    window.loadVcf_{jsid} = function(){{
      const btn = document.querySelector('#vload-{jsid} .load-btn');
      btn.disabled = true; btn.textContent = 'Loading…';
      document.getElementById('vstatus-{jsid}').style.display = '';
      setStatus(JSID, 'Downloading {jsid}.json…');

      const p1 = fetch('{jsid}.json')
        .then(r=>{{ if(!r.ok) throw new Error('HTTP '+r.status); return r.text(); }})
        .then(text=>{{
          setStatus(JSID,'Parsing variants (background thread)…');
          return parseJsonInWorker(text);
        }});

      const p2 = HAS_GFF
        ? fetch('{jsid}_features.json').then(r=>r.ok?r.json():[]).catch(()=>[])
        : Promise.resolve([]);

      Promise.all([p1, p2])
        .then(([data, features])=>{{
          document.getElementById('vload-{jsid}').style.display = 'none';

          _varData[JSID]  = data;
          _posIndex[JSID] = buildPosIndex(data);

          buildHeader(data, 'th-{jsid}', SNAMES, JSID);

          const dc = domChrom(data);
          const idx = _posIndex[JSID][dc] || [];
          const midEntry = idx[Math.floor(idx.length/2)];
          const midPos   = midEntry ? midEntry[0] : 0;
          if(midPos) {{
            renderWindow(JSID, dc, midPos, SNAMES);
          }} else {{
            renderAll(JSID, 'ALL', SNAMES);
          }}

          if(HAS_GFF && features.length){{
            document.getElementById('gff-wrap-{jsid}').style.display = '';
            setTimeout(()=>renderGffTrack(features, data, JSID), 50);
          }}
        }})
        .catch(err=>{{
          document.getElementById('vstatus-{jsid}').style.display = 'none';
          btn.disabled = false; btn.textContent = '↺ Retry';
          btn.insertAdjacentHTML('afterend',
            `<p class="muted" style="color:#f66;margin-top:.4rem">Error: ${{err.message}}</p>`);
        }});
    }};
  }})();
  </script>
</div>"""

def render_sample_page(sample_dir: Path, max_vcf_rows: int, fasta_name: str, gff_name: str,
                       index_url="index.html", generated_at="", has_gff=False):
    name = sample_dir.name
    fm, md_html = parse_md(sample_dir / "sample.md")
    sid  = re.sub(r"[^a-zA-Z0-9_-]", "_", name)

    meta_html = ""
    if fm:
        items = "".join(f"<dt>{_escape(k)}</dt><dd>{_escape(v)}</dd>" for k,v in fm.items())
        meta_html = f'<dl class="meta-dl">{items}</dl>'

    ov_html = _gallery_html(collect_images(sample_dir/"genome_overview"), "Genome overview", f"{sid}-ov")
    nv_html = _gallery_html(collect_images(sample_dir/"nucleotide_view"), "Nucleotide view",  f"{sid}-nv")

    vcfs = find_vcfs(sample_dir)
    vcf_metas = []
    vcf_blocks_html = ""
    if not vcfs:
        vcf_blocks_html = '<p class="muted">No VCF files found.</p>'
    else:
        for vcf_path in vcfs:
            jsid = re.sub(r"[^a-zA-Z0-9_]", "_",
                          f"{sid}__{vcf_path.stem.replace('.vcf','')}")
            vcf_blocks_html += _vcf_section_html(vcf_path, jsid, max_vcf_rows, fasta_name, gff_name, has_gff)
            _, col_headers, snames, csq_fmt, ann_fmt = read_vcf_header(vcf_path)
            if col_headers:
                vcf_metas.append((vcf_path, jsid, col_headers, snames, csq_fmt, ann_fmt))

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{_escape(name)} — Genomics Report</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;600;700&family=JetBrains+Mono:wght@400;700&display=swap" rel="stylesheet">
<script src="https://cdn.jsdelivr.net/npm/igv@2.15.5/dist/igv.min.js"></script>
<style>{{CSS}}</style>
</head><body>
<nav id="sidebar">
  <div class="sb-head"><div class="h1">Genomics Report</div>
  <div class="sub">{_escape(name)}</div>
  <div class="sub" style="margin-top:.15rem;font-size:10px;opacity:.6">{generated_at}</div></div>
  <div class="sb-search"><input type="text" placeholder="Filter samples…" oninput="sbFilter(this)"></div>
  <a class="sb-idx" href="{index_url}">⬆ Index</a>
  {{SIDEBAR_ITEMS}}
</nav>
<main id="main">
  <div style="margin-bottom:1rem"><a href="{index_url}" class="back-link">← Back to index</a></div>
  <h2 class="sample-name">{_escape(name)}</h2>
  {meta_html}
  {('<div class="sample-notes">'+md_html+'</div>') if md_html else ''}
  <div class="master-bar">
    <button class="master-btn" onclick="expandAll()">⊞ Expand all</button>
    <button class="master-btn" onclick="collapseAll()">⊟ Collapse all</button>
  </div>
  {ov_html}
  {nv_html}
  <div class="variants-section">
    <h3 class="variants-title">Variants
      <span class="count">({len(vcfs)} VCF{'s' if len(vcfs)!=1 else ''})</span>
    </h3>
    {vcf_blocks_html}
  </div>
</main>
<div id="lightbox">
  <button id="lb-close">✕</button>
  <img id="lb-img" src="" alt="">
  <div id="lb-cap"></div>
</div>
<div id="kbd-hint">Esc = close lightbox</div>
<div class="ts-footer">Generated: {generated_at}</div>
<script>{{JS}}</script>
</body></html>"""

    return html, vcf_metas

# ── Index page ────────────────────────────────────────────────────────────────

def render_index(samples, data_dir, generated_at=""):
    if not samples:
        return '<p class="muted">No sample directories found.</p>'
    colour_order, seen, sample_meta = [], set(), []
    for s in samples:
        fm,_ = parse_md(s/"sample.md")
        col = fm.get("colour", fm.get("color","")).strip().lower() or None
        if col not in seen: colour_order.append(col); seen.add(col)
        sample_meta.append((s, fm, col))
    colour_order.sort(key=lambda c: (c is None, c or ""))
    NAMED_BG = {"green":"rgba(39,174,96,.18)","red":"rgba(231,76,60,.18)",
                "yellow":"rgba(241,196,15,.15)","blue":"rgba(52,152,219,.18)",
                "orange":"rgba(230,126,34,.18)","purple":"rgba(155,89,182,.18)",
                "grey":"rgba(127,140,141,.18)","gray":"rgba(127,140,141,.18)"}
    rows_html = ""
    for col in colour_order:
        group = [(s,fm) for s,fm,c in sample_meta if c==col]
        if not group: continue
        if col:
            rows_html += (f'<tr><td colspan="6" style="padding:.3rem .75rem;font-size:11px;'
                          f'font-weight:700;text-transform:uppercase;letter-spacing:.06em;'
                          f'color:var(--muted);background:var(--surf2);'
                          f'border-bottom:1px solid var(--brd)">{col.capitalize()}</td></tr>')
        for s, fm in group:
            sid  = re.sub(r"[^a-zA-Z0-9_-]", "_", s.name)
            vcfs = find_vcfs(s)
            ov = (s/"genome_overview").is_dir() and any(p.suffix.lower() in IMAGE_EXTS for p in (s/"genome_overview").iterdir()) if (s/"genome_overview").is_dir() else False
            nv = (s/"nucleotide_view").is_dir()  and any(p.suffix.lower() in IMAGE_EXTS for p in (s/"nucleotide_view").iterdir())  if (s/"nucleotide_view").is_dir()  else False
            tick = lambda b: "✓" if b else "—"
            bg    = NAMED_BG.get(col, "rgba(128,128,128,.15)") if col else ""
            style = f"cursor:pointer;background:{bg}" if bg else "cursor:pointer"
            rows_html += (
                f'<tr onclick="window.location=\'{sid}.html\'" style="{style}">'
                f'<td><a href="{sid}.html" style="font-weight:600">{_escape(s.name)}</a></td>'
                f'<td>{_escape(fm.get("date",""))}</td>'
                f'<td>{_escape(fm.get("description",fm.get("desc","")))}</td>'
                f'<td class="tc">{tick(ov)}</td><td class="tc">{tick(nv)}</td>'
                f'<td class="tc">{len(vcfs) if vcfs else "—"}</td>'
                f'</tr>')
    ts = (f'<p class="muted" style="margin-top:.75rem;font-size:11px">Generated: {generated_at}</p>'
          if generated_at else "")
    return f"""<div class="idx-section">
  <h2>Sample index <span class="count">({len(samples)} samples)</span></h2>
  <table class="idx-table">
  <thead><tr><th>Sample</th><th>Date</th><th>Description</th>
  <th class="tc">Overview</th><th class="tc">Nucleotide</th><th class="tc">VCFs</th></tr></thead>
  <tbody>{rows_html}</tbody></table>{ts}</div>"""

# ── Assembly ──────────────────────────────────────────────────────────────────

def _sidebar_items(samples, current_sid=None):
    out = []
    for s in samples:
        sid = re.sub(r"[^a-zA-Z0-9_-]", "_", s.name)
        cls = "sb-item" + (" current" if sid==current_sid else "")
        out.append(f'<a class="{cls}" href="{sid}.html">{_escape(s.name)}</a>')
    return "\n".join(out)

def _wrap(raw, samples, current_sid):
    return (raw.replace("{CSS}", CSS)
               .replace("{JS}",  JS)
               .replace("{SIDEBAR_ITEMS}", _sidebar_items(samples, current_sid)))

def _render_worker(args):
    sample_dir, max_vcf_rows, fasta_name, gff_name, generated_at, has_gff = args
    html, vcf_metas = render_sample_page(sample_dir, max_vcf_rows, fasta_name, gff_name,
                                          generated_at=generated_at, has_gff=has_gff)
    return sample_dir, html, vcf_metas

def build_report(samples, data_dir, out_dir, max_vcf_rows, gff_path=None, fasta_name=None, skip_json=False, ncpu=1):
    out_dir.mkdir(parents=True, exist_ok=True)
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    has_gff = gff_path is not None and gff_path.exists()
    gff_name = gff_path.name if has_gff else None
    tasks = [(s, max_vcf_rows, fasta_name, gff_name, generated_at, has_gff) for s in samples]

    if ncpu > 1 and len(samples) > 1:
        print(f"Rendering HTML ({len(samples)} samples, {ncpu} workers)…", file=sys.stderr)
        with Pool(processes=ncpu) as pool:
            results = list(tqdm(pool.imap(_render_worker, tasks),
                                total=len(tasks), desc="HTML", unit="sample", file=sys.stderr))
    else:
        results = list(tqdm((_render_worker(t) for t in tasks),
                            total=len(tasks), desc="HTML", unit="sample", file=sys.stderr))

    for sample_dir, html, _ in results:
        sid = re.sub(r"[^a-zA-Z0-9_-]", "_", sample_dir.name)
        (out_dir / f"{sid}.html").write_text(_wrap(html, samples, sid), encoding="utf-8")

        sample_out = out_dir / sample_dir.name
        sample_out.mkdir(exist_ok=True)
        for view_dir in ["nucleotide_view", "genome_overview"]:
            src = sample_dir / view_dir
            if src.is_dir():
                dst = sample_out / view_dir
                if not dst.exists():
                    try: dst.symlink_to(src.resolve())
                    except Exception: pass
    
    all_metas = [(sd, m) for sd, _, metas in results for m in metas]
    if all_metas and not skip_json:
        print(f"\nStreaming {len(all_metas)} VCF(s) → JSON…", file=sys.stderr)
        for sample_dir, (vcf_path, jsid, col_headers, snames, csq_fmt, ann_fmt) in \
                tqdm(all_metas, desc="VCF→JSON", unit="vcf", file=sys.stderr):
            json_path = out_dir / f"{jsid}.json"
            try:
                n, _ = write_variant_json(vcf_path, json_path, col_headers, snames,
                                       csq_fmt, ann_fmt, max_vcf_rows)
                print(f"  {vcf_path.name}: {n} variants → {json_path.name}", file=sys.stderr)
            except Exception as e:
                print(f"  [WARN] {vcf_path.name}: {e}", file=sys.stderr)
                continue

            if has_gff:
                coord_ranges = get_vcf_coord_range(json_path)
                all_features = []
                for chrom, (start, end) in coord_ranges.items():
                    pad = max(5000, (end - start) // 10)
                    feats = extract_gff_features(gff_path, chrom, max(1, start-pad), end+pad)
                    all_features.extend(feats)
                feat_path = out_dir / f"{jsid}_features.json"
                feat_path.write_text(json.dumps(all_features, separators=(',',':')),
                                     encoding="utf-8")
    elif skip_json:
        print("\n[INFO] --skip-json active: Bypassing variant JSON conversion and GFF compilation.", file=sys.stderr)

    idx_body = render_index(samples, data_dir, generated_at)
    sb       = _sidebar_items(samples)

    changelog_html = ""
    changelog_path = data_dir / "index.md"
    if changelog_path.exists():
        _, md_html = parse_md(changelog_path)
        if md_html.strip():
            changelog_html = (
                f'<div class="section" style="margin-bottom:1.5rem">'
                f'<div class="section-hdr" onclick="toggleSection(\'index-notes\')">'
                f'<h3>Notes</h3>'
                f'<span class="section-arrow open" id="ar-index-notes">▶</span></div>'
                f'<div class="section-body open" id="sb-index-notes" '
                f'style="padding:1rem 1.25rem;font-size:13px;line-height:1.7">'
                f'{md_html}</div></div>'
            )

    idx_page = f"""<!DOCTYPE html><html lang="en"><head>
<meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Genomics Report — {_escape(data_dir.name)}</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;600;700&family=JetBrains+Mono:wght@400;700&display=swap" rel="stylesheet">
<style>{CSS}</style></head><body>
<nav id="sidebar">
  <div class="sb-head"><div class="h1">Genomics Report</div>
  <div class="sub">{len(samples)} sample{"s" if len(samples)!=1 else ""}</div>
  <div class="sub" style="margin-top:.15rem;font-size:10px;opacity:.6">{generated_at}</div></div>
  <div class="sb-search"><input type="text" placeholder="Filter samples…" oninput="sbFilter(this)"></div>
  <a class="sb-idx" href="index.html">⬆ Index</a>{sb}
</nav>
<main id="main">{changelog_html}{idx_body}</main>
<div id="lightbox"><button id="lb-close">✕</button><img id="lb-img" src="" alt=""><div id="lb-cap"></div></div>
<script>{JS}</script></body></html>"""
    (out_dir / "index.html").write_text(idx_page, encoding="utf-8")

    total_mb = sum(f.stat().st_size for f in out_dir.iterdir() if f.is_file()) / 1e6
    print(f"\nDone → {out_dir}/  ({total_mb:.1f} MB total)", file=sys.stderr)


def main():
    args     = parse_args()
    data_dir = Path(args.data)
    out_dir  = Path(args.output_dir)
    gff_path = Path(args.gff) if args.gff else None
    fasta_path = Path(args.fasta) if args.fasta else None
    
    if not data_dir.exists():
        print(f"[ERROR] '{data_dir}' not found.", file=sys.stderr); sys.exit(1)
    if gff_path and not gff_path.exists():
        print(f"[WARN] GFF not found: {gff_path} — continuing without GFF track.", file=sys.stderr)
        gff_path = None
    if fasta_path:
        if not fasta_path.exists():
            print(f"[WARN] FASTA not found: {fasta_path}. IGV track may fail.", file=sys.stderr)
            fasta_path = None
        else:
            out_dir.mkdir(parents=True, exist_ok=True)
            link_target = out_dir / fasta_path.name
            if not link_target.exists():
                try:
                    link_target.symlink_to(fasta_path.resolve())
                    fai_path = Path(str(fasta_path) + ".fai")
                    if fai_path.exists():
                        (out_dir / fai_path.name).symlink_to(fai_path.resolve())
                except Exception as e:
                    print(f"[WARN] Failed to symlink FASTA: {e}", file=sys.stderr)

    if gff_path and gff_path.exists():
        out_dir.mkdir(parents=True, exist_ok=True)
        link_target = out_dir / gff_path.name
        if not link_target.exists():
            try:
                link_target.symlink_to(gff_path.resolve())
                # Expose tabix index to IGV if it exists
                for ext in [".tbi", ".csi"]:
                    idx_path = Path(str(gff_path) + ext)
                    if not idx_path.exists():
                        idx_path = gff_path.with_suffix(gff_path.suffix + ext)
                    if idx_path.exists():
                        (out_dir / idx_path.name).symlink_to(idx_path.resolve())
            except Exception as e:
                print(f"[WARN] Failed to symlink GFF: {e}", file=sys.stderr)

    samples = discover_samples(data_dir)
    for sample_id in samples:
        # Assuming BAM files live in data_dir/sample_id/sample_id.bam
        sample_dir = data_dir / sample_id
        bam_candidate = sample_dir / f"{sample_id}.bam"
        
        if bam_candidate.exists():
            # Create a dedicated web-accessible data directory inside report/
            web_sample_dir = out_dir / "data" / sample_id
            web_sample_dir.mkdir(parents=True, exist_ok=True)
            
            target_bam = web_sample_dir / f"{sample_id}.bam"
            target_bai = web_sample_dir / f"{sample_id}.bam.bai"
            
            try:
                if not target_bam.exists():
                    target_bam.symlink_to(bam_candidate.resolve())
                
                # Check for standard .bam.bai or .bai variants
                source_bai = Path(str(bam_candidate) + ".bai")
                if not source_bai.exists():
                    source_bai = bam_candidate.with_suffix(".bai")
                    
                if source_bai.exists() and not target_bai.exists():
                    target_bai.symlink_to(source_bai.resolve())
            except Exception as e:
                print(f"[WARN] Failed to symlink alignment tracks for {sample_id}: {e}", file=sys.stderr)    
    if not samples:
        print("[WARN] No sample directories found.", file=sys.stderr)
    print(f"Found {len(samples)} sample(s) in '{data_dir}'", file=sys.stderr)
    
    build_report(samples, data_dir, out_dir, args.max_vcf_rows,
                 gff_path=gff_path, fasta_name=fasta_path.name if fasta_path else None, skip_json=args.skip_json, ncpu=args.ncpu)

if __name__ == "__main__":
    main()