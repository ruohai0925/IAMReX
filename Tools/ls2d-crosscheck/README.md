# ls2d-crosscheck — material from cross-validating IAMReX against ls2d-baseline (2026-08-25)

Moved here from `taskpreparation/ls2d-baseline/` so that ls2d-baseline contains only the teaching code and its
comparison with the papers. VTK/npy snapshots in `out_archive/` are git-ignored; each run directory keeps its `diag.csv`, `summary.txt` and log.

Before running the Python scripts, activate a virtual environment and install
[`requirements.txt`](requirements.txt) from this directory:

```bash
python3 -m pip install -r requirements.txt
```

| item | content |
|---|---|
| `PLAN.md` | the original plan / milestone log of ls2d-baseline (Chinese) |
| `deviations.md` | D1–D14: every place where IAMReX's level-set branch differs from Sussman et al. 1999 / Sussman–Fatemi 1999, with status. Led to IAMReX PR #77 (D4 `mass_fix` delta support on the wrong field, D7 docs, D13 Taylor–Green pressure comment) and PR #78 (D14 `ns.epsilon_massfix` + phase-volume diagnostic) |
| `validation.md` | full chronological validation log incl. code-to-code tables ls2d vs IAMReX (Taylor–Green, Rayleigh–Taylor, RSV area drift), multigrid fix history and ablation studies |
| `algorithm_with_iamrex_map.md` | equation ↔ ls2d ↔ IAMReX function correspondence table |
| `compare_iamrex.py`, `read_amrex_plotfile.py` | read AMReX plotfiles (Header / Level_0 / Cell_D) into numpy; overlay interface contours from IAMReX and ls2d VTK |
| `out_archive/` | all ls2d run directories from the development phase: RSV cross-checks (`rsv4_*`), IAMReX-parameter RT (`rt_iamrex64`, `rt_code_to_code.png`), multigrid tuning (`rtq_*`, `tg64_*`), reinit ablations (`vtx_*`, `zal_*`, `disk_*`), breaking-wave trial (`bw128`) |

IAMReX runs used for the cross-checks live in `Tutorials/RSV/run_*` (RSV, `mass.txt` columns 5–6 = phase volume) and were
made with `Tutorials/RSV/amr2d.gnu.ex` (USE_MPI=FALSE, `prob.probtype` switch for Taylor–Green / RT).
