# Benchmark Tools

This directory has a small potential benchmark, GP fixture generation, and a
couple of sweep scripts for comparing CPU, GPU scalar, and GPU batched paths.

## Build

Build the potential benchmark explicitly:

```bash
cmake -S . -B build
cmake --build build --target potential_benchmark.e
```

For GPU runs, configure the build with the GPU backend:

```bash
cmake -S . -B build-sycl -D GPU_BACKEND=sycl
cmake --build build-sycl --target potential_benchmark.e
```

The binary lands in the build directory, for example
`build/potential_benchmark.e` or `build-sycl/potential_benchmark.e`.

## Potential Benchmark

Show options:

```bash
./build/potential_benchmark.e --benchmark-help
```

List available interaction and external potentials:

```bash
./build/potential_benchmark.e --benchmark-list
```

Example interaction benchmark:

```bash
./build/potential_benchmark.e \
  --benchmark-kind interaction \
  --benchmark-potential aziz \
  --benchmark-method all \
  --benchmark-iterations 1000000
```

Example external benchmark:

```bash
./build/potential_benchmark.e \
  --benchmark-kind external \
  --benchmark-potential free \
  --benchmark-method value
```

Non-benchmark flags are passed through to the normal PIMC setup parser. Defaults
are supplied for the usual required simulation options, so short benchmark
commands work without a full input deck.

Common options:

- `--benchmark-method value|gradient|pair|action|all`
- `--benchmark-samples N`
- `--benchmark-min R`
- `--benchmark-max R`
- `--benchmark-check-batch`
- `--benchmark-gpu`
- `--benchmark-check-gpu`

Use `--benchmark-gpu` or `--benchmark-check-gpu` only with a GPU build.

## GP Fixtures

Generate synthetic GP input data:

```bash
python3 tools/benchmarks/generate_gp_fake_data.py \
  --output-dir gp_bench_data/points-4096 \
  --points 4096
```

This writes `gp_input.ini` and `gp_training.dat`. These files are just benchmark
fixtures.

## Sweep Scripts

`run_gp_action_sweep.py` uses `potential_benchmark.e` to time GP action calls.
It builds the benchmark, creates missing GP fixtures, runs the requested cases,
and writes:

- `benchmark_results/gp_action_sweep.csv`
- `benchmark_results/gp_action_sweep.md`

Example:

```bash
python3 tools/benchmarks/run_gp_action_sweep.py \
  --gp-points 1024,4096 \
  --particles 16,64 \
  --slices 16,64 \
  --iterations 16
```

`run_gp_pimc_sweep.py` runs end-to-end PIMC cases with the GP potential. It
builds `pimc.e`, creates missing GP fixtures, and writes:

- `benchmark_results/gp_pimc_sweep.csv`
- `benchmark_results/gp_pimc_sweep.md`

Example:

```bash
python3 tools/benchmarks/run_gp_pimc_sweep.py \
  --gp-points 16384 \
  --particles 16,64 \
  --slices 16,64 \
  --updates default,bisection
```

Both scripts default to running CPU, GPU scalar, and GPU batched modes. Useful
flags:

- `--cpu-build-dir PATH`
- `--gpu-build-dir PATH`
- `--data-root PATH`
- `--output-dir PATH`
- `--compiler PATH`
- `--boost-root PATH`
- `--gpu-selector SELECTOR`
- `--skip-build`
- `--skip-cpu`
- `--skip-gpu-scalar`
- `--skip-gpu-batched`
- `--dry-run`

The scripts set `PIMC_ACTION_CALL_STATS=1` and
`PIMC_BATCHED_POTENTIAL_STATS=1`. GPU runs use `ONEAPI_DEVICE_SELECTOR` from
`--gpu-selector`; GPU scalar mode sets `PIMC_DISABLE_BATCHED_POTENTIAL=1`.
