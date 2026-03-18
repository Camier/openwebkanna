### 2026-03-14 / Codex
- hypothesis: with `MolDetv2` held constant as the detector, backend yield differs materially between `MolGrapher` and `MolScribe` on the same extraction set.
- change: add a small comparative OCSR harness and run it on `Capps`, `Stafford`, and `Gericke`.
- metrics: per-backend detector crops, validated SMILES, invalid SMILES, unique valid SMILES, and per-paper overlap.
- success criteria: produce a repeatable artifact that makes backend tradeoffs explicit without relying on terminal-only evidence.
- result: `MolGrapher` produced `41` unique valid SMILES vs `40` for `MolScribe` on the same `48` MolDetv2 crops, with lower invalid output (`5` vs `8`). `MolScribe` won `Stafford` (`4` valid vs `3`), but `MolGrapher` stayed ahead on `Capps` and `Gericke`. Cross-backend overlap stayed low (`jaccard_avg=0.111`), so they are not interchangeable.
- decision: keep both backends available, keep `MolGrapher` as the safer default on this panel, and treat the harness as the baseline before any routing or backend-default change. Publish the exact backend disagreement sets in `disagreements.json` and summarize a few examples in `report.md` so routing work can target specific chemistry mismatches instead of aggregate counts alone.
- follow-up: `run-chemical-ocsr-disagreement-analysis.sh` now buckets those mismatches with RDKit-backed categories. On the current panel, the dominant pattern is `substituent_drift` (`15`) followed by `coverage_gap` (`11`), with only a small number of `placeholder_or_attachment` (`2`) and `radical_or_query_atom` (`1`) cases.

### 2026-03-14 / Codex
- hypothesis: backend complementarity is large enough that a conservative fusion view is more actionable than a single default backend, but only if the merged pool is tiered by agreement risk.
- change: add a fusion analysis pass on top of `run.json` + `disagreement_analysis.json` and classify candidates into `consensus`, `complementary`, `review_required`, and `risky`.
- metrics: union upper bound, auto-keep count (`consensus + complementary`), review-required count, risky count, and per-paper tier split.
- success criteria: produce a repeatable artifact that quantifies both the upside of backend fusion and the portion that still needs human review.
- result: the RDKit-valid union reaches `71` candidates, but only `21` are conservative auto-keep (`10` consensus + `11` complementary coverage gains). `46` candidates are structurally competing variants that need review, and `4` are risky placeholder/radical-style artifacts. `Gericke` drives most review load (`42` review-required of `54` union candidates), while `Capps` is mostly pure coverage gain (`8` auto-keep of `11` union candidates).
- decision: keep `MolGrapher` as the default extractor, but treat backend fusion as the review queue generator. The new `fusion_candidates.json` / `fusion_report.md` artifacts make that split explicit, and `review_queue.json` plus `review_queue/*.md` now expose the per-paper operator queue directly instead of pretending the backend union can be ingested blindly.
- follow-up: add flat adjudication exports (`review_queue.csv`, `review_queue.jsonl`, and per-paper `review_queue/*.csv` / `review_queue/*.jsonl`) so the queue can be handed to human review or a backend triage worker without scraping Markdown.
- next step: add a post-adjudication materializer that reads the queue back in, preserves explicit human decisions when present, and emits `accepted` / `rejected` / `pending` exports for downstream ingestion or curation workflows.
- implementation note: keep manual-review and defaulted-materialization outputs separate (`adjudication_manual_*` vs `adjudication_defaulted_*`) so conservative auto-materialization never overwrites the explicit-review lane.
- follow-up: emit a compact `accepted_catalog` alongside the full accepted exports so downstream ingestion can consume `effective_smiles` with page/block/crop provenance without carrying the full review queue schema.
