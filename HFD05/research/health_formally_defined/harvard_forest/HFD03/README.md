# HFD03 — Reconstruction calibration, failure localization and discriminative validation

Additive revision; HFD01/HFD02S frozen at4053fa. [Closure](reports/Closure.md), [contract](reports/hfd03_contract.json), [protocol](config/protocol.json), [referee map](reports/PaperV_Revision_Map.md), [figures](figures/README.md).

Use Python3.9+ with NumPy2.0.2, preferably the existing private HFD02S runtime. No Lean rebuild, new worktree or copied cache required.

```sh
python scripts/verify_hfd03_harvard_forest.py
python scripts/verify_hfd03_harvard_forest.py --cache
python scripts/verify_hfd03_harvard_forest.py --replay
```

Portable audit reads committed evidence. Cache audit reuses read-only HFD01 inputs and HFD03 isolated identity/sanitized caches. Full replay regenerates only HFD03 outputs under the frozen design; no original HFD02S bank is written.

Private generated masked CSVs and compressed missingness sufficient banks are in .runs/de7f3568b1942f9b; public extra identity inputs in .inputs. These are not copied historical scientific banks. Source IDs, namespace/provenance rules, seed algorithm, all candidate IDs, all mask assignments and sufficient MC counts are committed. The manifests bind implementation and output hashes.

Scientific limits are first-class results: coarse DBH undercoverage; entry ambiguity; incomparable exposed E2 anchors; model-dependent discrimination; no mass-near-threshold synthetic cases; unresolved broad-grid cells. Operational Harvard Forest Health remains NOT_AUTHORIZED.
