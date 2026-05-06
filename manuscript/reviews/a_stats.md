## Verdict: FIX-THEN-SHIP

## Critical (must fix before submission)

- **§3.5, sensitivity-sweep, random-negatives line.** Manuscript claims "V2-316M − V1 paired Δ = +0.069 [+0.062, +0.088] under random negatives". Reproduction from `results/probes/layer2__Geneformer-V2-316M__neg-random.parquet` and `layer2__geneformer__neg-random.parquet` gives Δ = +0.0149 [+0.003, +0.027] (n=375, seed=0, 2000 boots; identical at 5000 boots). The +0.069 number is wrong by a factor of ~4.6×. This is the only sentence in the sensitivity section that doesn't reproduce, and it's a load-bearing one (it argues the V1→V2 ordering survives random-negative sampling). The qualitative claim "V1→V2 ordering held" is still true with the correct +0.015 (CI excludes zero), but the magnitude is much smaller than under matched negatives, which actually changes the narrative slightly (random negatives narrow, not amplify, the V1→V2 gap, because the easier task compresses the dynamic range). Fix the number, and consider rewording "Crucially the V1→V2 ordering and the within-V2 null both held" to acknowledge the gap shrinks under easier negatives.

- **§3.1, L3 trained-vs-random-init CIs (V2-104M and V2-316M).** Reproduction with seed=0, 2000 bootstrap resamples on per-fold paired differences:
  - V2-104M − V2-104M-rand: Δ = +0.179, CI [+0.143, +0.215]. Manuscript: "Δ = +0.18 [+0.12, +0.26]". Both endpoints disagree (lower 0.14 vs 0.12; upper 0.22 vs 0.26).
  - V2-316M − V2-316M-rand: Δ = +0.256, CI [+0.214, +0.297]. Manuscript: "Δ = +0.26 [+0.23, +0.33]". Both endpoints disagree (lower 0.21 vs 0.23; upper 0.30 vs 0.33).
  These look like they were computed by subtracting per-model independent CIs rather than bootstrapping the paired difference. Same point estimate, wrong CI. Either (a) recompute using paired bootstrap to match the rest of the paper, or (b) caveat that these are unpaired CI ranges.

- **§3.1, L3 V2-316M − V1 CI.** Manuscript: "Δ = -0.020 [-0.044, +0.026], p = 0.625". Reproduction: Δ = -0.020, CI [-0.053, +0.013], p = 0.625. Point estimate and p reproduce; CI does not (manuscript is narrower and shifted right by ~0.013). A reviewer recomputing from the parquets will flag.

- **Methods §2 vs implementation: bootstrap resample count.** Methods says "5,000 resamples". `morpheus/analysis/figures_v2.py:181` uses 2000. The number is small enough that 5000 vs 2000 doesn't move CIs noticeably (verified), but a methods reviewer who runs the code will note the discrepancy. Either change Methods to "2,000 resamples" or change the implementation to match.

## Nice-to-have (would strengthen but not blocking)

- **§3.1, V2-316M − bag CI upper bound.** Manuscript [+0.13, +0.29]; reproduction [+0.13, +0.27] (mean Δ = +0.199). The +0.29 in the abstract ("CI [+0.13, +0.29] across all three scales") is correctly the V1−bag bound, but writing the same [+0.13, +0.29] for V316−bag specifically in §3.1 line 480 misstates the upper bound by 0.02.

- **§3.3 layer-wise wording.** "layer 9 was the strongest individual layer at 0.430, and the final-layer readout we use elsewhere tied for best" — the parquet has layer 9 = 0.4296 and layer 12 = 0.4301. Layer 12 is the unique (rounded-to-three-digits-tied) best, not layer 9. Either say "layer 12 was the strongest at 0.430, with layer 9 tied (0.430)" or just "layer 12 was best (0.430)." Current phrasing reads as if layer 9 wins outright.

- **§3.4 RSA random-init range.** Manuscript: "Random-init transformers were similar (-0.001 to -0.003)". Actual values: V1-rand +0.0010, V104-rand -0.0022, V316-rand -0.0026. Range is more like [-0.003, +0.001]. Trivially small numerically but the stated sign of the upper end is wrong.

- **Abstract L2 p-value rounding.** Abstract says "p=10⁻¹⁶"; actual p = 2.4×10⁻¹⁶. Convention varies but writing "p<10⁻¹⁵" or "p=2×10⁻¹⁶" is more honest about the order.

- **§3.5 mean vs median labelling.** L2 MLP comparisons explicitly use "median AUC" (0.587, 0.701). L3 MLP comparison "V2-316M L3 MLP 0.567 (down from 0.698)" silently switches to means (median is 0.571; mean is 0.567). Add "(mean)" to keep it consistent.

- **Methods §2 selectivity claim.** "median 0.500 at every scale, mean within ±0.01 of 0.5" — verified: V1 perm mean 0.495, V104 perm 0.490, V316 perm 0.496. Mean for V2-104M (0.490) is technically just outside ±0.01 of 0.5 (it's 0.0105 below). Tighten to "±0.015" or recheck.

## Verified clean

- §3.2 L2 paired contrasts: V1→V2-104M Δ=+0.054 [+0.041,+0.068], p=2.4×10⁻¹⁶; V1→V2-316M Δ=+0.051 [+0.038,+0.066], p=3.2×10⁻¹⁵; V2-316M − V2-104M Δ=−0.003 [−0.014,+0.008], p=0.45; V2-316M − bag Δ=+0.167 [+0.151,+0.182], p=4×10⁻⁵⁰. All reproduce exactly.
- §3.3 L1 V2-316M − bag Δ=+0.019 [−0.001,+0.043], p=0.31. Reproduces.
- §3.3 L1 V2-316M − V2-104M Δ=−0.007 [−0.034,+0.017], p=1.00. Reproduces.
- §3.3 L1 means: V1 0.354, V104 0.430, V316 0.424, bag 0.404. All match `stats_summary.json`.
- §3.3 L1 random-init means: V1-rand 0.219, V104-rand 0.191, V316-rand 0.149. Match.
- §3.1 L3 V2-316M − V2-104M Δ=−0.011 [−0.028,+0.007], p=0.44. Reproduces (matches earlier reframe fix).
- §3.1 L3 V1−bag CI [+0.13,+0.29], V104−bag [+0.14,+0.27]. Reproduce.
- §3.1 V1-rand mean 0.42, V104-rand mean 0.53, V316-rand mean 0.44 (L3). Match.
- §3.4 RSA: V1 ρ=−0.028, V104 ρ=−0.004, V316 ρ=+0.022, bag ρ=+0.026. All match parquet exactly.
- §3.5 sensitivity DoRothEA-only V316−V1 Δ=+0.067 [+0.054,+0.084] and CollecTRI V316−V1 Δ=+0.056 [+0.041,+0.073]. Reproduce.
- §3.5 MLP medians: V1 L2 0.587, V316 L2 0.701, V316 L3 mean 0.567, all L1 MLP means. Reproduce.
- §3.1 sign-test: walked back from (1/2)⁶=1/64 to (1/2)³=1/8 on the three independent trained-vs-random-init contrasts only; explicit non-combination caveat present for the trained-vs-bag contrasts citing shared baseline + shared trained-AUC components. Honest and correct.
- n=5 floor: explicitly noted in §2 ("smallest two-sided Wilcoxon p-value attainable is 0.0625, which we report unmodified") and §4.2. No fold-level test reports p<0.0625. ✓
- "Fails to exceed" / "fails to reject" language: consistently paired with caveats about not having performed equivalence testing (§3.3, §4.2). No "matches"/"is equivalent to" slip-ups in the final tex.
- Paired Wilcoxon labelling: every "paired" claim in §3 is on per-fold or per-TF differences; no Mann-Whitney conflation found.
- Cell-type count 39 (was 47), parameter ratio 31× consistent across §3.1 and §4.1, "did not exceed bag-of-genes" wording — all confirmed in the rendered tex.
