# TMLR submission variant

Anonymized double-blind submission package for *Transactions on Machine Learning Research*.

## What changed vs the bioRxiv variant

- Document class switched from `article` to `\usepackage{tmlr}` (official TMLR style file from `JmlrOrg/tmlr-style-file`).
- `\author{}` block: replaced with `Anonymous Authors / Anonymous Affiliation`. The `tmlr` package suppresses the author block in submission mode anyway, but the placeholder makes the intent explicit and survives a stray switch to `[preprint]` mode.
- All identifying URLs (GitHub repository, Zenodo concept and version DOIs, ORCID, contact email, institutional affiliation, fancyhdr running banner) removed from the body.
- Acknowledgments, Funding, and AI-assistance paragraph rewritten to avoid identifying information; vendor name removed from the AI-assistance disclosure. The substance of the disclosure (LLM-assisted preparation; author responsibility for content) is preserved.
- `\bibliographystyle{plainnat}` → `\bibliographystyle{tmlr}`. Self-citations: none — there is no prior work by the author to anonymize in the bibliography.
- Figure paths localized to `figures/` so the variant compiles standalone.

## What did NOT change

- Body text (Introduction, Methods, Results, Discussion).
- Numerical results, figures, and statistical claims.
- Bibliography (`references.bib`) — same 42 entries.
- Supplementary Table S1 contents — only the title-block author info was redacted.

## Compiling

```bash
tectonic manuscript.tex
tectonic supplementary.tex
```

## Switching to camera-ready or preprint mode

Submission (default, anonymous):
```latex
\usepackage{tmlr}
```
Camera-ready (accepted, author identities visible, "Published in TMLR" banner):
```latex
\usepackage[accepted]{tmlr}
```
Non-anonymous preprint posting (author identities visible, no TMLR banner):
```latex
\usepackage[preprint]{tmlr}
```

After accepted-state switch, also restore: title-block author, the github/Zenodo URLs in `Reproducibility` and `Code availability`, the original `Acknowledgments` paragraph, the original `Funding` statement, and the vendor name in the AI-assistance disclosure if desired.

## Bundling for OpenReview

OpenReview accepts a single anonymized PDF for the main manuscript, plus one supplementary file (PDF or ZIP, up to 100 MB). For this submission:
- **Main:** `manuscript.pdf`
- **Supplementary:** `supplementary.pdf` (Table S1 with per-comparison statistics)

A code/data ZIP can be added as a second supplementary file if desired; see the TMLR author guide.
