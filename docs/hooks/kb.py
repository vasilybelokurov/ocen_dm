"""MkDocs hook for the oCen_dm knowledge base.

* Serves files that live outside ``docs/``: ``plots/*.png`` at ``plots/``,
  JOURNAL.md as month pages under ``journal/``, OCEN_DM_TAILS_PROJECT.md as
  ``spec.md``.
* Generates two pages from committed sources so they cannot drift:
  ``kb/data/sources.md`` from ``provenance/datasets.yaml`` and
  ``kb/literature.md`` from ``docs/kb/literature.json``.
* Rewrites links in every page against the file's real location in the repo:
  targets inside ``docs/`` or ``plots/`` stay site links; other tracked repo
  files (code, configs, tests) point to GitHub; untracked targets
  (``results/``, ``data/``) become plain text, since collaborators do not have
  them. Absolute local links into the repository (Codex reviews) are mapped
  the same way.
"""
from __future__ import annotations

import json
import os
import posixpath
import re
import subprocess
from collections import defaultdict
from pathlib import Path

import yaml
from mkdocs.structure.files import File

ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT/"docs"
GITHUB = "https://github.com/vasilybelokurov/ocen_dm/blob/main/"
SCIX_ID = "RSXq3-NtTTCApvfPUbdwrg"
LINK = re.compile(r"(!?)\[((?:[^\[\]]|\[[^\]]*\])*)\]\((<[^>]+>|[^)\s]+)(\s+\"[^\"]*\")?\)")
SPACED = re.compile(r"(!?)\[([^\[\]]*)\]\(([^()<>\"]*\s[^()<>\"]*)\)")
_REAL = {}          # virtual src_uri -> real path of a generated page (for link resolution)
_TRACKED = None


def _tracked():
    global _TRACKED
    if _TRACKED is None:
        out = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True).stdout
        _TRACKED = set(out.split("\n"))
    return _TRACKED


# ---------------------------------------------------------------- generated pages
def _journal_pages():
    text = (ROOT/"JOURNAL.md").read_text()
    head, *entries = re.split(r"(?m)^(?=## )", text)
    months = defaultdict(list)
    for e in entries:
        m = re.match(r"## (\d{4}-\d{2})", e)
        months[m.group(1) if m else "undated"].append(e)
    pages = {}
    for month, es in months.items():
        pages[f"journal/{month}.md"] = f"# Journal, {month}\n\n" + "".join(es)
    index = ["# Journal", "",
             "Chronological record of the project, split by month from the repository's `JOURNAL.md`.",
             "Later entries can supersede earlier ones; corrections are recorded as new entries.", ""]
    index += [f"- [{m}]({m}.md): {len(months[m])} entries" for m in sorted(months, reverse=True)]
    pages["journal/index.md"] = "\n".join(index)+"\n"
    return pages


def _sources_page():
    reg = yaml.safe_load((ROOT/"provenance/datasets.yaml").read_text())["datasets"]
    rows = ["# Registered data sources", "",
            "Generated at build time from `provenance/datasets.yaml`, the checked registry of external data.",
            "`verified` means the link was checked live and recorded in the registry; unverified entries say why.",
            "Inputs not yet registered are listed on [Known gaps](gaps.md).", "",
            "| Key | Source | Role | Paper | Data | Files | Verified |", "|---|---|---|---|---|---|---|"]
    for key, d in reg.items():
        paper = d.get("paper") or {}
        data = d.get("data") or {}
        links = [f"[arXiv]({paper['arxiv']})" if paper.get("arxiv") else "",
                 (f"[journal]({paper['journal']})" if str(paper.get("journal", "")).startswith("http")
                  else str(paper.get("journal", ""))),
                 f"[DOI](https://doi.org/{paper['doi']})" if paper.get("doi") else "",
                 f"[ADS]({paper['ads']})" if paper.get("ads") else ""]
        dlink = (f"[{data['doi']}](https://doi.org/{data['doi']})" if data.get("doi") else
                 f"[{data.get('kind', 'link')}]({data['record_url']})" if data.get("record_url") else
                 f"[{data.get('kind', 'link')}]({data['url']})" if data.get("url") else data.get("kind", ""))
        files = ", ".join(f"`{f['name']}`" for f in d.get("files") or [] if isinstance(f, dict) and f.get("name"))
        ver = "yes" if paper.get("verified", True) and data.get("verified", True) else "no"
        role = " ".join(str(d.get("role", "")).split())
        rows.append(f"| `{key}` | {d.get('authors', '')}: {d.get('title', '')} | {role} | "
                    f"{' '.join(x for x in links if x)} | {dlink} | {files} | {ver} |")
    return "\n".join(rows)+"\n"


def _literature_page():
    d = json.loads((DOCS/"kb/literature.json").read_text())
    papers = d["papers"]
    topics = ["data-source", "omega-cen-mass-imbh", "cluster-dynamics-methods", "tides-and-dm-survival",
              "tails-and-streams", "mw-potential-bar-lmc", "stream-fitting-methods", "stellar-populations-mf", "other"]
    names = {"data-source": "Data sources", "omega-cen-mass-imbh": "Omega Cen mass, remnants and central mass",
             "cluster-dynamics-methods": "Cluster dynamics methods (DF, Jeans, multimass, N-body)",
             "tides-and-dm-survival": "Tides and dark-matter survival", "tails-and-streams": "Tails and streams",
             "mw-potential-bar-lmc": "Milky Way potential, bar and LMC", "stream-fitting-methods": "Stream-fitting methods",
             "stellar-populations-mf": "Stellar populations and mass function", "other": "Other"}
    n_lib = sum(1 for p in papers if p.get("bibcode") and p["relevance"] != "tangential")
    out = ["# Literature", "",
           f"!!! abstract \"SciX library\"\n    All {n_lib} papers below with a bibcode and relevance *core* or *context* are in the public library "
           f"[omega Cen dark matter (oCen_dm)](https://scixplorer.org/public-libraries/{SCIX_ID}) "
           f"([ADS view](https://ui.adsabs.harvard.edu/public-libraries/{SCIX_ID})).", "",
           "Generated from `docs/kb/literature.json`, compiled from every citation in the repository "
           f"({d.get('generated', '')}); identifiers were checked in ADS. *Core*: the project uses or depends on it; "
           "*context*: cited for background; *tangential*: mentioned once (not in the library).", ""]
    for t in topics:
        sel = sorted((p for p in papers if p["topic"] == t), key=lambda p: ({"core": 0, "context": 1}.get(p["relevance"], 2), p["cite"]))
        if not sel:
            continue
        out += [f"## {names[t]}", "", "| Paper | Title | Relevance | Links |", "|---|---|---|---|"]
        for p in sel:
            links = []
            if p.get("bibcode"):
                links.append(f"[ADS](https://ui.adsabs.harvard.edu/abs/{p['bibcode'].replace('&', '%26')})")
            if p.get("arxiv"):
                links.append(f"[arXiv](https://arxiv.org/abs/{p['arxiv']})")
            if p.get("doi"):
                links.append(f"[DOI](https://doi.org/{p['doi']})")
            if not links:
                links.append("unresolved")
            title = (p.get("title") or "").replace("|", "/")
            out.append(f"| {p['cite']} | {title} | {p['relevance']} | {' '.join(links)} |")
        out.append("")
    bad = d.get("erroneous_identifiers_exclude") or []
    if bad:
        out += ["## Identifiers found to be wrong", "",
                "These identifiers appear in repository text (mostly Codex research reports) but point to other papers; "
                "the correct identifiers are used above.", "", "| As written | Note |", "|---|---|"]
        for b in bad:
            ident = b.get("identifier") or b.get("arxiv") or b.get("doi") or b.get("bibcode") or ""
            out.append(f"| `{ident}` | {' '.join(str(b.get('note') or b.get('reason') or '').split())} |")
    return "\n".join(out)+"\n"


def on_files(files, config):
    for rel in sorted(t for t in _tracked() if t.startswith("plots/") and t.endswith(".png")):   # committed plots only
        files.append(File(rel, str(ROOT), config["site_dir"], config["use_directory_urls"]))
    pages = dict(_journal_pages())
    pages["spec.md"] = (ROOT/"OCEN_DM_TAILS_PROJECT.md").read_text()
    pages["kb/data/sources.md"] = _sources_page()
    pages["kb/literature.md"] = _literature_page()
    for uri, content in pages.items():
        _REAL[uri] = ROOT/("JOURNAL.md" if uri.startswith("journal/") else
                           "OCEN_DM_TAILS_PROJECT.md" if uri == "spec.md" else "docs/"+uri)
        files.append(File.generated(config, uri, content=content))
    return files


# ---------------------------------------------------------------- link rewriting
def _rewrite(page_uri, target):
    """Return a new link target, or None to turn the link into plain text."""
    raw = target[1:-1] if target.startswith("<") else target
    if re.match(r"^[a-z][a-z0-9+.-]*:", raw, re.I) or raw.startswith("#"):
        return target
    path, frag = (raw.split("#", 1)+[""])[:2]
    line = re.search(r":(\d+)(?:-\d+)?$", path)
    if line and not Path(path).suffix.startswith(f".{line.group(1)}"):
        path, frag = path[:line.start()], frag or f"L{line.group(1)}"
    real_page = _REAL.get(page_uri, DOCS/page_uri)
    if path.startswith("/"):
        p = Path(path)
        try:
            real = ROOT/p.relative_to(ROOT)
        except ValueError:
            return None
    elif not path:
        return target
    else:
        real = Path(os.path.normpath(real_page.parent/path))
        if not real.exists() and (ROOT/path).exists():      # repo-root-relative links in reviews
            real = ROOT/path
    try:
        rel_repo = real.relative_to(ROOT).as_posix()
    except ValueError:
        return None
    if rel_repo.startswith("docs/") and real.exists() and not rel_repo.startswith("docs/hooks/"):
        virtual = rel_repo[len("docs/"):]
    elif rel_repo.startswith("plots/") and rel_repo in _tracked():
        virtual = rel_repo
    elif rel_repo == "JOURNAL.md":
        virtual = "journal/index.md"
    elif rel_repo == "OCEN_DM_TAILS_PROJECT.md":
        virtual = "spec.md"
    elif rel_repo in _tracked() or (real.is_dir() and any(t.startswith(rel_repo+"/") for t in _tracked())):
        kind = "tree/main/" if real.is_dir() else ""
        url = (GITHUB.replace("blob/main/", kind) if kind else GITHUB)+rel_repo
        return url+(f"#{frag}" if frag else "")
    else:
        return None
    new = posixpath.relpath(virtual, posixpath.dirname(page_uri) or ".")
    if re.fullmatch(r"L\d+(-L?\d+)?", frag):      # source-line anchors only exist on GitHub
        frag = ""
    return new+(f"#{frag}" if frag else "")


def on_page_markdown(markdown, page, config, files):
    uri = page.file.src_uri

    def sub(m):
        bang, text, target, title = m.group(1), m.group(2), m.group(3), m.group(4) or ""
        new = _rewrite(uri, target)
        if new is None:
            shown = target.strip("<>")
            return f"*[figure not in the repository: `{shown}`]*" if bang else (text if text.strip() else f"`{shown}`")
        return f"{bang}[{text}]({new}{title})"

    def spaced(m):          # e.g. "[packet:2021](question packet)": not a path, keep the text
        return m.group(2)

    # leave fenced code blocks untouched
    parts = re.split(r"(?ms)(^```.*?^```)", markdown)
    return "".join(p if p.startswith("```") else LINK.sub(sub, SPACED.sub(spaced, p)) for p in parts)
