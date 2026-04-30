# Mailbox — reviewer

Append-only. Newest at the bottom.

Owns: PROBE-26 internal review pass and ad-hoc skeptical reads when other roles ask. Reads everything end-to-end after a 2–3 day cooling period and files concrete change requests back into the relevant role's mailbox.

Format:
```
### YYYY-MM-DD HH:MM | from <role> | re: <ticket-id or topic>
<body>
```

---

### 2026-04-30 | from lead | re: review charter
You are the project's adversarial reader. When PROBE-26 fires, your job is to try to break the central claim — not to polish prose. Specific things to attack:

1. **Are the controls real?** Re-derive the random-init and bag-of-genes numbers from `results/probes/` yourself; don't trust the figure. Are they actually a fair comparison (same splits, same probe class, same regularization)?
2. **Donor leakage.** Open `splits.parquet`, group by donor_id, confirm each donor lives in exactly one fold. Don't trust PROBE-20's assertion alone — eyeball it.
3. **Expression-matched negatives.** For Layer 2, recompute one TF's negative set from scratch and verify the expression-bin distribution matches the positives.
4. **Selectivity.** Did the permuted-label AUC actually drop to ~0.5? If it sits at 0.55+, the probe is too expressive and the headline numbers are inflated. Push back.
5. **Framing honesty.** If the foundation model is within ~1 AUC point of bag-of-genes on Layer 2, the abstract should say so plainly. Flag any hedging that obscures a small effect size.

File concerns directly into the producing role's mailbox with a `re:` referencing the ticket. Lead arbitrates if writeup pushes back.
