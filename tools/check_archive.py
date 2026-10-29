"""Check the English archive package; this is not a theorem prover."""

from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

from pypdf import PdfReader


ROOT = Path(__file__).resolve().parents[1]
PAPERS = tuple(path.name for path in sorted((ROOT / "papers").iterdir()) if (path / "main.tex").is_file())
SKIP = {".git", "build", "tmp", "__pycache__"}
TEXT_SUFFIXES = {".md", ".tex", ".bib", ".wls", ".ps1", ".py", ".json", ".txt"}
checks = []


def require(name, condition, detail=""):
    checks.append({"check": name, "passed": bool(condition), "detail": detail})


paper_records = []
for slug in PAPERS:
    folder = ROOT / "papers" / slug
    source = "\n".join(path.read_text(encoding="utf-8-sig") for path in sorted(folder.glob("*.tex")))
    bib = (folder / "references.bib").read_text(encoding="utf-8-sig")
    require(f"{slug}: blank TeX byline", re.search(r"\\author\s*\{\s*\}", source))
    keys = re.findall(r"@\w+\s*\{\s*([^,]+),", bib)
    duplicates = [key for key, count in Counter(keys).items() if count > 1]
    require(f"{slug}: unique bibliography keys", not duplicates, str(duplicates))
    cited = set()
    for group in re.findall(r"\\cite\w*\*?(?:\[[^]]*\])*\{([^}]+)\}", source):
        cited.update(key.strip() for key in group.split(","))
    missing = sorted(cited - set(keys))
    require(f"{slug}: bibliography coverage", not missing, str(missing))
    labels = re.findall(r"\\label\{([^}]+)\}", source)
    require(f"{slug}: unique labels", len(labels) == len(set(labels)))
    reader = PdfReader(folder / "manuscript.pdf")
    text = " ".join((page.extract_text() or "") for page in reader.pages)
    text = re.sub(r"\s+", " ", text)
    author = str((reader.metadata or {}).get("/Author", "")).strip()
    require(f"{slug}: blank PDF author", not author, repr(author))
    # PDF text extraction can omit visually present interword spacing.
    compact_text = re.sub(r"\s+", "", text)
    for phrase in ("I genuinely do not know", "Fable 5", "DeepSeek", "GPT", "human expert review"):
        require(f"{slug}: disclosure '{phrase}'", phrase.replace(" ", "").casefold() in compact_text.casefold())
    require(f"{slug}: no unresolved PDF references", "??" not in text)
    require(f"{slug}: no placeholder byline", "name(s) to be supplied" not in text)
    require(f"{slug}: English PDF", re.search(r"[\u3400-\u9fff]", text) is None)
    log = (folder / "build" / "main.log").read_text(encoding="utf-8", errors="replace")
    errors = re.findall(r"^!.*|.*(?:undefined references|Citation .* undefined|Reference .* undefined|Overfull \\[hv]box).*", log, re.MULTILINE)
    require(f"{slug}: clean final LaTeX log", not errors, "\n".join(errors))
    paper_records.append({"paper": slug, "pages": len(reader.pages), "pdf_author": author, "citation_count": len(cited)})

files = []
for path in sorted(ROOT.rglob("*")):
    rel = path.relative_to(ROOT)
    if not path.is_file() or any(part in SKIP for part in rel.parts) or rel.as_posix() == "archive-checks.json":
        continue
    raw = path.read_bytes()
    files.append({"path": rel.as_posix(), "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()})
    if path.suffix in TEXT_SUFFIXES or path.name in {".gitignore", ".gitattributes"}:
        text = raw.decode("utf-8-sig")
        require(f"{rel.as_posix()}: English text", re.search(r"[\u3400-\u9fff]", text) is None)
        private_path = re.search(r"[A-Za-z]:[/\\]Users[/\\]", text)
        require(f"{rel.as_posix()}: no local user paths", private_path is None)
        token = re.search(r"gh[pousr]_[A-Za-z0-9_]{20,}|-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----", text)
        require(f"{rel.as_posix()}: no common credential patterns", token is None)

report = {
    "checked_at_utc": datetime.now(timezone.utc).isoformat(),
    "kind": "English manuscript build and archive checks",
    "human_expert_review": False,
    "historical_Wolfram_scripts_rerun": False,
    "papers": paper_records,
    "checks": checks,
    "files": files,
}
(ROOT / "archive-checks.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
failed = [check for check in checks if not check["passed"]]
print(json.dumps({"checks": len(checks), "failed": failed, "papers": paper_records}, indent=2))
raise SystemExit(1 if failed else 0)
