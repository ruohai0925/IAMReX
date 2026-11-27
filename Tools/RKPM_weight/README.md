# RKPM Weight Tool

## Architecture

See [`ARCHITECTURE.md`](ARCHITECTURE.md) for the maintained English code
architecture diagrams, module responsibilities and runtime data contracts. The
architecture document must be reviewed and updated in the same change whenever
the code structure, runtime flow, interfaces or MPMD protocol changes.

`main.py` separates the weight solver backend from the output transport through
two independent runtime arguments:

| Argument | Values | Purpose |
|---|---|---|
| `--solver` | `rkpm` / `ml` | Traditional RKPM solve or Transolver inference |
| `--transport` | `file` / `mpmd` | Write `.id/.lag` files or return weights directly to IAMReX |
| `--batch-size` | positive integer | Markers processed per vectorized batch (default: `4096`) |
| `--device` | `cpu` / `cuda` / `auto` | PyTorch device used by the ML backend (default: `cpu`) |

Both solver backends use batched array operations. Traditional RKPM assembles
all polynomial bases and moment matrices in a batch and solves the independent
10-by-10 systems with batched NumPy linear algebra. The ML backend performs
coordinate preprocessing and Transolver inference in batches. In MPMD mode,
the fixed stencil metadata is packed and validated once, then reused on every
exchange; weights remain arrays throughout the solve and send path.

`batch_size` limits temporary memory without changing marker or stencil order.
The default is intended as a balanced starting point. Reduce it when the Python
rank has limited memory; increase it only after benchmarking the target CPU or
accelerator. Numerical results are independent of how markers are split into
batches, up to normal floating-point roundoff.

The default ML device remains `cpu` for portability and reproducibility. Use
`device = cuda` on a GPU allocation to require CUDA, or `device = auto` to use
CUDA when PyTorch can see a GPU and otherwise fall back to CPU. `--device` has
no effect on the traditional NumPy RKPM backend.

## Python Environment

Run the tool with a Python 3 interpreter in which the required dependencies are
installed. In the commands below, `python3` means that interpreter; activate a
virtual environment first or replace `python3` with its absolute path when
needed.

From `Tools/RKPM_weight/`, install the standard dependencies from
[`requirements.txt`](requirements.txt):

```bash
python3 -m pip install -r requirements.txt
```

NumPy is required for all modes. Point-cloud generation and visualization also
require Matplotlib, ML inference requires PyTorch, and MPMD transport requires
`mpi4py`. The MPI implementation used by `mpi4py` must be compatible with the
one used to build IAMReX.

## Ellipsoid Surface Markers

`generate_ellipsoid_points.py` generates a requested number of surface markers
using E3 repulsion (NumPy and SciPy required). Both the command line and Python
function interpret the three semi-axis lengths directly as x/y/z, in the
supplied order. To put the long axis along x, supply its length first, as in
`--axes 0.0655 0.0437 0.0437`.

```bash
python3 generate_ellipsoid_points.py -n 512 \
  --axes 0.0655 0.0437 0.0437 --output ellipsoid_points.txt --plot-areas
```

For a transverse ellipsoid with its long axis along y, use
`--axes 0.0437 0.0655 0.0437 --mirror-axis y`. This requires an even marker
count and reflects the first half across the plane `y = center_y` in the
same order. E3 relaxation includes the full mirrored cloud in its energy.
Generate the `.id` and `.lag` from the new coordinate file.

The coordinate file uses `# X Y Z` followed by exactly N coordinate rows, ready
for `main.py --geometry ... --body-frame`. The center defaults to the origin.
Adjacent `.areas.csv` and `.areas.json` files report estimated surface Voronoi
areas and their deviation from S/N; `--plot-areas` also writes a PNG and requires
Matplotlib. All reports use the same final x/y/z semi-axes as the coordinates.
Match `particle_inputs.radius`, `radius2`, and `radius3` to those final body-frame
semi-axes when configuring an ellipsoid.

## Runtime Configuration File

All Python-side runtime arguments can be stored in an AMReX-style `key = value`
configuration file and loaded with one option. The bundled
`unit_test/inputs.rkpm` is a test-case example; production runs should use a
case-specific configuration file:

```bash
python3 main.py --config /path/to/case/inputs.rkpm
```

Paths in the configuration are resolved relative to the configuration file,
not the shell's working directory. Explicit command-line arguments override
the corresponding configuration values, for example:

```bash
python3 main.py --config /path/to/case/inputs.rkpm \
  --solver ml --device cuda --batch-size 2048 --check-action abort
```

The configuration controls the Python weight service. The IAMReX executable,
MPI rank counts and `IAMReX_MPMD=1` remain part of the `mpiexec` command.

## File Mode

Build new stencils and mappings from a point cloud:

```bash
python3 main.py \
  --solver rkpm --transport file \
  --inputs /path/to/case/inputs.3d \
  --geometry /path/to/case/lagrangian_points.txt \
  --body-frame \
  --output-prefix /path/to/case/rkpm_mappings
```

Point-cloud generation requires an explicit coordinate frame. Choose exactly
one of the following options:

- Use `--body-frame` when the point-cloud coordinates are relative to the body.
  The tool reads `particle_inputs.x`, `particle_inputs.y` and
  `particle_inputs.z` from the AMReX inputs file and translates the points to
  that world-coordinate center before constructing their grid stencils.
- Use `--world-frame` when the point-cloud coordinates are already expressed in
  the AMReX world coordinates bounded by `prob_lo` and `prob_hi`. No translation
  is applied.

Omitting both options is an error; the tool does not guess the coordinate
frame. In a runtime configuration file, `body_frame = true` selects body-frame
coordinates and `body_frame = false` explicitly selects world-frame
coordinates.

For ellipsoid point clouds, pass the semi-axes of the marker surface using
`--axes A B C` or `axes = A B C` in the runtime configuration file. Use the
same length units as the point cloud and grid. The tool divides the estimated
ellipsoid surface area equally among the markers. Omitting `--axes` estimates
the semi-axes from half the input cloud's span in each direction before applying
`--angle`. This estimate assumes the ellipsoid's principal axes align with the
input coordinate axes and the cloud samples near the axis endpoints. Supply
the semi-axes explicitly if the input cloud is already tilted.

Read marker positions and stencils from existing `.id/.lag` files and recompute
only the ML weights:

```bash
python3 main.py \
  --solver ml --transport file \
  --inputs /path/to/case/inputs.3d \
  --id /path/to/case/rkpm_mappings.id \
  --lag /path/to/case/rkpm_mappings.lag \
  --model-dir /path/to/model_directory \
  --model-code /path/to/transolver_slim.py \
  --output-prefix /path/to/case/rkpm_mappings_ml
```

The traditional RKPM backend can use the same mapping-recompute path by changing
`--solver ml` to `--solver rkpm`.

When building a new mapping from point-cloud geometry, all 3x3x3 support cells
are constructed by broadcasting integer grid offsets. The generator therefore
scales with the number of markers and their 27-point supports instead of
repeatedly scanning the full Eulerian bounding grid for every marker.

## MPMD Mode

C++ still reads `rkpm_mappings.id/.lag` during initialization. The `.lag` file
provides the fixed `i/j/k/Vcell/eps` stencil. The Python server receives marker
positions from C++, computes 27 weights per marker and returns only the weights:

```bash
mpiexec -n 2 env IAMReX_MPMD=1 \
    /path/to/amr3d.gnu.MPI.ex /path/to/case/inputs.3d \
  : -n 1 python3 \
      /path/to/IAMReX/Tools/RKPM_weight/main.py \
      --config /path/to/case/inputs.rkpm
```

The example uses two CFD ranks and one Python-server rank; adjust only the CFD
rank count for the target case because the Python server currently requires
exactly one rank. Set `solver` in the case configuration, or override it with
`--solver rkpm` or `--solver ml`, to select the backend without changing the
MPMD protocol.

Marker positions use double precision (`MPI_DOUBLE`) so cell-index validation
does not lose coordinate precision in transit. Returned weights remain single
precision (`MPI_FLOAT`) to match the ML model output. Coordinates within
floating-point roundoff of an internal grid line are snapped to that line
before applying `floor`, so the exact line consistently belongs to the cell on
its high-index side.

Before sending MPMD weights, the server can check the zeroth moment and the
first moment in cell units. `check_action` selects `off`, `warn`, or `abort`,
and `check_interval` controls how often the check runs. The check uses the
actual float32 array sent to C++. With the provided model and mapping, the
traditional RKPM backend satisfies the default `1e-6` tolerances, while the ML
backend currently exceeds the first-moment tolerance and therefore reports a
warning (or aborts when `check_action = abort`).

## Current MPMD Protocol Limitation

The protocol returns weights but not stencil `i/j/k` indices. A marker may move
inside the center cell of its original stencil, but it cannot cross a cell
boundary. The server aborts the MPI job if it detects such a crossing, preventing
C++ from applying new weights to stale cell indices. Supporting continuously
moving particles requires extending both the C++ and Python sides to transfer
updated stencil indices.

## Tests

Run the complete unit-test suite from any working directory with:

```bash
/path/to/IAMReX/Tools/RKPM_weight/unit_test/run_all_tests.sh
```

See `unit_test/README.md` for dependency handling, interpreter selection and
details of the individual tests.

The tool root intentionally contains only `src/`, `unit_test/`, `main.py`,
`README.md` and `ARCHITECTURE.md`. Bundled mapping/input fixtures, the standalone
alignment verifier and generated reference images live under `unit_test/`.
