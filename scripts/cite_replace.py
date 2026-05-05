#!/usr/bin/env python3
"""Convert inline (Author Year) citations in manuscript.tex to \\citep{key} / \\citet{key}.

Reads /Users/williamcole/Documents/Development/chemistry/Morpheus/manuscript/manuscript.tex
Writes the same file, replacing inline citation strings.
"""
import re
from pathlib import Path

SRC = Path("/Users/williamcole/Documents/Development/chemistry/Morpheus/manuscript/manuscript.tex")

# Citation key map. Order matters — more specific patterns first.
# Each entry: (regex matching the inline citation, citation key).
# The regex must capture two groups when textual: name, year. For parenthetical, just match the whole inline.

# Use a simple registry: pattern → key. Then we do two passes:
# Pass A — parenthetical citations: \(([A-Za-z][^)]*?\d{4}[a-z]?)\)
#   → look up the matched body against our patterns, replace with \citep{key}
# Pass B — textual citations: Author( and Author| et al\.)? \(YYYY\)
#   → \citet{key}

# Map: list of (sub-regex on inline body, key)
KEY_MAP = [
    (r"theodoris.*?2023", "theodoris2023"),
    (r"cui.*?2024", "cui2024"),
    (r"hao.*?2024", "hao2024scfoundation"),
    (r"schaar.*?2025", "schaar2025nicheformer"),
    (r"rosen.*?2023", "rosen2023uce"),
    (r"boiarsky.*?2025", "boiarsky2025"),
    (r"kedzierska.*?2025", "kedzierska2025"),
    (r"ahlmann-eltze[^()]*2025", "ahlmanneltze2025perturbation"),
    (r"ahlmann-eltze.*?2023", "ahlmanneltze2023normalization"),
    (r"ahlmann.*?huber.*?2023", "ahlmanneltze2023normalization"),
    (r"csendes.*?2025", "csendes2025"),
    (r"helical.*?2024", "helical2024"),
    (r"venkatesh.*?2026", "venkatesh2026"),
    (r"qiu.*?2025", "qiu2025bioLLM"),
    (r"civale.*?2026", "civale2026"),
    (r"garcia-alonso.*?2019", "garciaalonso2019dorothea"),
    (r"m.{0,4}ller-?dott.*?2023", "mullerdott2023collectri"),
    (r"müller-?dott.*?2023", "mullerdott2023collectri"),
    (r"badia-i-mompel.*?2022", "badiaimompel2022decoupler"),
    (r"czi cell science program.*?2025", "cellxgene2025"),
    (r"cellxgene.*?2025", "cellxgene2025"),
    (r"wolf.*?2018", "wolf2018scanpy"),
    (r"wolf.*?2020", "wolf2020transformers"),
    (r"pedregosa.*?2011", "pedregosa2011sklearn"),
    (r"paszke.*?2019", "paszke2019pytorch"),
    (r"vaswani.*?2017", "vaswani2017attention"),
    (r"devlin.*?2019", "devlin2019bert"),
    (r"bommasani.*?2021", "bommasani2021foundation"),
    (r"diehl.*?2016", "diehl2016cellontology"),
    (r"osumi-sutherland.*?2021", "osumisutherland2021hca"),
    (r"lun.*?2016", "lun2016scran"),
    (r"conneau.*?2018", "conneau2018"),
    (r"hewitt.*?2019", "hewitt2019"),
    (r"hewitt.*?liang.*?2019", "hewitt2019"),
    (r"belinkov.*?2022", "belinkov2022"),
    (r"kriegeskorte.*?2008", "kriegeskorte2008rsa"),
    (r"zhang.*?bowman.*?2018", "zhang2018random"),
    (r"pratapa.*?2020", "pratapa2020beeline"),
    (r"reimers.*?gurevych.*?2019", "reimers2019sbert"),
    (r"reimers.*?2019", "reimers2019sbert"),
    (r"wilcoxon.*?1945", "wilcoxon1945"),
    (r"efron.*?1979", "efron1979bootstrap"),
    (r"zimmerman.*?2021", "zimmerman2021pseudoreplication"),
    (r"abdelaal.*?2019", "abdelaal2019"),
    (r"luecken.*?2022", "luecken2022scib"),
    (r"kaplan.*?2020", "kaplan2020scaling"),
    (r"hoffmann.*?2022", "hoffmann2022chinchilla"),
    (r"chinchilla", "hoffmann2022chinchilla"),
    (r"alain.*?bengio.*?2016", "alainbengio2016probes"),
    # Second-author-only fallback patterns (for textual cites like "Anders 2025")
    (r"^anders.*?2025", "ahlmanneltze2025perturbation"),
    (r"^liang.*?2019", "hewitt2019"),
    (r"^huber.*?2023", "ahlmanneltze2023normalization"),
    (r"^langefeld.*?2021", "zimmerman2021pseudoreplication"),
    (r"^bengio.*?2016", "alainbengio2016probes"),
    (r"^bowman.*?2018", "zhang2018random"),
    (r"^liang.*?2017", "hewitt2019"),
    (r"^toutanova.*?2019", "devlin2019bert"),
    (r"^teichmann.*?2021", "osumisutherland2021hca"),
]


def lookup_key(body: str) -> str | None:
    body_lower = body.lower()
    for pat, key in KEY_MAP:
        if re.search(pat, body_lower):
            return key
    return None


def replace_parenthetical(text: str) -> tuple[str, int, list[str]]:
    """Replace `(Author Year)` with `\\citep{key}` where the body matches a known author/year."""
    # Match parens containing a year, e.g. (Theodoris et al., 2023) or (Hewitt and Liang, 2019)
    pattern = re.compile(r"\(([^()]+\b\d{4}[a-z]?\b[^()]*)\)")
    n = 0
    misses = []

    def sub(m: re.Match) -> str:
        nonlocal n
        body = m.group(1)
        key = lookup_key(body)
        if key:
            n += 1
            return f"\\citep{{{key}}}"
        misses.append(body)
        return m.group(0)

    return pattern.sub(sub, text), n, misses


def replace_textual(text: str) -> tuple[str, int, list[str]]:
    """Replace `Author et al. (YYYY)` and `Author and Author (YYYY)` with `\\citet{key}`."""
    # Several patterns; do them in order.
    n = 0
    misses = []

    # Pattern 1: Author1 and Author2 (YYYY)  — e.g. "Hewitt and Liang (2019)"
    p1 = re.compile(r"\b([A-Z][A-Za-z\-]+ and [A-Z][A-Za-z\-]+) \((\d{4}[a-z]?)\)")
    # Pattern 2: Author et al. (YYYY)
    p2 = re.compile(r"\b([A-Z][A-Za-z\-]+ et al\.?) \((\d{4}[a-z]?)\)")
    # Pattern 3: Single Author (YYYY) — careful, can over-match
    p3 = re.compile(r"\b([A-Z][A-Za-z\-]+) \((\d{4}[a-z]?)\)")

    for p in (p1, p2, p3):
        def sub(m: re.Match) -> str:
            nonlocal n
            body = f"{m.group(1)} {m.group(2)}"
            key = lookup_key(body)
            if key:
                n += 1
                return f"\\citet{{{key}}}"
            misses.append(body)
            return m.group(0)
        text = p.sub(sub, text)

    return text, n, misses


def main() -> None:
    text = SRC.read_text()
    # Pre-process: collapse linebreaks inside potential citation parens (Author\nYYYY).
    # This is a common pandoc-induced wrapping that breaks the regex matcher.
    text = re.sub(r"\(([^()]*?)\n([^()]*?\b\d{4}[^()]*?)\)", lambda m: f"({m.group(1).strip()} {m.group(2).strip()})", text)
    text, n_paren, misses_p = replace_parenthetical(text)
    text, n_textual, misses_t = replace_textual(text)
    SRC.write_text(text)
    print(f"replaced {n_paren} parenthetical citations and {n_textual} textual citations")
    if misses_p or misses_t:
        print("UNMATCHED CITATIONS (need manual review):")
        for m in misses_p:
            print(f"  paren: ({m})")
        for m in misses_t:
            print(f"  textual: {m}")


if __name__ == "__main__":
    main()
