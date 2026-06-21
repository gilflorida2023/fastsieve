# fastsieve — Configurable Wheel Segmented Prime Sieve

**fastsieve** counts and enumerates primes in large ranges using a
segmented Sieve of Eratosthenes with configurable wheel factorization
(mod 2, 6, 30, 210, or 2310) and 128-bit integers.  It can sieve up to
~10³⁸, limited only by time.

## Algorithm

### Sieve of Eratosthenes (Segmented)

The classic sieve marks multiples of each known prime as composite.
Rather than allocating a bitmap for the entire range (impossible for
large N), **segmentation** processes the number line in fixed-size
windows.  Only primes up to √N need to be stored in memory — about
50 million primes for N = 10¹⁸, which fits in ~400 MB.

For each window:
1. Initialize the sieve buffer (all candidates = prime)
2. For each stored base prime, walk its multiples within the window
   and mark them composite
3. Scan the window for survivors — each survivor is a new prime
4. Buffer discovered primes in memory; flush to state file in 65K-entry batches

### Configurable Wheel Factorization

Numbers divisible by the first k primes can never be prime (except
those primes themselves).  The wheel pre-excludes them: of every
W consecutive integers (where W is the primorial), only φ(W) residues
are coprime to W and need to be checked.

| Wheel (--wheel) | Primorial | Primes in Wheel | Residues | % Candidates | % Removed |
|---|---|---|---|---|---|
| 2 | 2 | {2} | 1 | 50% | 50% |
| 6 | 2×3 | {2,3} | 2 | 33.3% | 66.7% |
| 30 | 2×3×5 | {2,3,5} | 8 | 26.7% | 73.3% |
| **210** (default) | **2×3×5×7** | **{2,3,5,7}** | **48** | **22.9%** | **77.1%** |
| 2310 | 2×3×5×7×11 | {2,3,5,7,11} | 480 | 20.8% | 79.2% |

The default `--wheel 210` offers the best balance of candidate reduction
vs. overhead on modern CPUs.  Larger wheels (2310) reduce candidates
further but increase modulus overhead and LUT size.

The wheel residues for the default mod-210 wheel are:

```
1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
101, 103, 107, 109, 113, 121, 127, 131,
137, 139, 143, 149, 151, 157, 163, 167, 169, 173,
179, 181, 187, 191, 193, 197, 199, 209
```

These are all numbers < 210 not divisible by 2, 3, 5, or 7.

### Dynamic Base Primes

Rather than pre-sieving all primes up to √N at startup, the base
prime list grows via **bootstrapping by squaring**:

1. Start with `{2, 3, 5, 7}`
2. Sieve the range `[8, 49]` using these — now we have primes up to 49
3. Sieve `[50, 2401]` using primes up to 49 — now we have primes up to 2401
4. Sieve `[2402, 5,764,801]` — primes up to 5.7M
5. Continue until we have primes up to √(current segment end)

Each level squares the maximum known prime, so 5 levels suffice for
any 128-bit target.

### 128-bit Integers

Standard C `uint64_t` overflows when computing `p * p` for
p > 2³² (~4.3×10⁹), limiting the safe target to ~10¹⁹.  **fastsieve**
uses `unsigned __int128` (GCC/Clang extension), which allows
`p * p` up to p = 2⁶⁴−1 and targets up to ~3.4×10³⁸.

### Origin of the Wheel Sieve

- **Pritchard (1981)**: "A Sublinear Additive Sieve for Finding Prime
  Numbers" — first formal description of wheel factorization.
- **Tomás Oliveira e Silva (2001)**: Bucket sieve algorithm for
  cache-efficient segmentation at extreme ranges.
- **Kim Walisch (primesieve)**: State-of-the-art implementation using
  a mod-210 wheel, SIMD, and bucket sieving; achieves 0.4 seconds
  for π(10⁹) on a modern CPU.

## Compilation

```sh
gcc -O3 -march=native -o fastsieve fastsieve.c
```

Requires GCC or Clang on a 64-bit platform (for `__int128` support).

## Usage

```sh
./fastsieve [--wheel N] [buffer_size] [target] [-c] [-s] [-r] [-R] [-o file]
```

| Argument | Description |
|---|---|---|
| `--wheel N` | Wheel modulus: 2, 6, 30, **210** (default), 2310 |
| `buffer_size` | Segment window in natural numbers (optional, auto-optimized to ~143K) |
| `target` | Sieve up to this number (required for sieve, report, and resume) |
| `-c` | Count only — no state file, faster for large targets |
| `-s` | Save state file — overrides `-c`; uses batched writes (minimal I/O penalty) |
| `-r` | Report from existing `primes_state.bin` (no sieving) |
| `-R` | Resume from checkpoint — load `primes_state.ckpt`, continue sieving from the last committed position (requires target > last sieved). State file is appended to automatically (unless `-c`). |
| `-o file` | Write all discovered primes to a file (use `-` for stdout) |

### Examples

| Command | What it does |
|---|---|
| `./fastsieve 100000000` | Sieve to 100M (target), save state file |
| `./fastsieve -c 1000000000000` | Count primes ≤ 10¹², no disk writes |
| `./fastsieve -c -o primes.txt 1000000000` | Count + write all primes to file |
| `./fastsieve -c -s 1000000000000` | Count + save state file (batched writes) |
| `./fastsieve -c -s -o primes.txt 1000000000` | State file + text output in count mode |
| `./fastsieve -r` | Print summary from existing state file |
| `./fastsieve -r -o primes.txt` | Reconstruct prime list from state file |
| `./fastsieve -s 10000000` | Sieve to 10M, create state + checkpoint files |
| `./fastsieve -R 20000000` | Resume from 10M → 20M, append to state file |
| `./fastsieve -R -c 50000000` | Resume, count only — no state writes, fastest |
| `./fastsieve -R 5000000` | Target ≤ last sieved → no-op, prints existing count |
| `./fastsieve -R 50000000` *(after crash)* | Resume from last committed checkpoint; re-processes ≤10 segments |
| `./fastsieve --wheel 30 1000000000 -c` | Count primes with mod-30 wheel |
| `./fastsieve --wheel 2310 1000000000 -c` | Count primes with mod-2310 wheel |
| `./fastsieve -c -r` | Error (mutually exclusive) |
| `./fastsieve -R -r` | Error (mutually exclusive) |
| `./fastsieve` | Print help |

## Performance

Benchmarked on an **AMD Ryzen 5 PRO 8500GE** (Zen 4c, 6 cores/12 threads,
single-threaded run, GCC 14.2.0 `-O3 -march=native`):

| Target | π(target) | Time | Segments |
|---|---|---|---|
| 10⁶ | 78,498 | 0.004 sec | 7 |
| 10⁷ | 664,579 | 0.03 sec | 70 |
| 10⁸ | 5,761,455 | 0.35 sec | 699 |
| 10⁹ | 50,847,534 | 3.9 sec | 6,983 |
| 10¹⁰ | 455,052,511 | 48 sec | 69,823 |
| 10¹¹ | ~4.1B | ~8 min* | 698,255 |
| 10¹² | ~37.6B | ~1.4 hr* | 6,982,555 |

\* Estimated from 10¹⁰ scaling (algorithm is O(N log log N)).

The inner marking loop uses incremental residue/block tracking instead of
128-bit modulo and division, giving a **~2.3× speedup** over a naive
implementation.  State file writes use batched I/O (65K entries per flush),
so the performance difference between default and `-c` is small.

## Resume & Crash Recovery

When running with state saving (`-s` or default mode), fastsieve maintains
a **checkpoint file** (`primes_state.ckpt`) that records the last committed
position.  This enables resuming after an interruption.

### How it works

1. **`primes_state.bin`** — append-only entry list.  Never overwritten, always
   consistent.  Entries are flushed in batches (65K per write) for performance.
2. **`primes_state.ckpt`** — 48-byte checkpoint with `last_sieved`, `total_count`,
   `entry_count`, `original_target`, and a magic number.  Updated atomically
   after every 10 segments and at final completion.

### Atomic update

The checkpoint is written to a temporary file first, then moved into place with
POSIX `rename()`, which is atomic on the same filesystem.  If the process is
killed during the write, the previous checkpoint survives — no corruption.

### What you lose on a crash

At most **10 segments** (≈1.4M numbers with the default buffer).  The next
resume starts from the last committed checkpoint and re-processes those
segments.  The state file may contain garbage bytes beyond `entry_count`;
these are overwritten on resume.

### Completing a partial run

```sh
# Example: start a long sieve
./fastsieve -s 1000000000000

# … machine crashes or Ctrl+C …
# Later, resume from wherever it left off:
./fastsieve -R 1000000000000
```

The `original_target` in the checkpoint tells you what the first run aimed
for.  If you accidentally specify a smaller target, fastsieve prints a warning
but continues.  If the new target is ≤ the last sieved position, it's a no-op
(the count is printed).

## State File Format

`primes_state.bin` stores discovered primes for resumability.
Each entry is 32 bytes:

```
offset 0: prime_lo (uint64_t, little-endian)
offset 8: prime_hi (uint64_t, little-endian)
offset 16: next_lo (uint64_t, little-endian)
offset 24: next_hi (uint64_t, little-endian)
```

Where `prime = (hi << 64) | lo` and `next` is the first multiple of
`prime` that falls in or after the current segment.

### Checkpoint File

`primes_state.ckpt` (64 bytes) stores the last committed resume position.
Updated atomically every 10 segments.

```
offset  0: last_sieved_lo (uint64_t, little-endian)
offset  8: last_sieved_hi (uint64_t, little-endian)
offset 16: total_count (uint64_t)
offset 24: entry_count (uint64_t)
offset 32: original_target_lo (uint64_t, little-endian)
offset 40: original_target_hi (uint64_t, little-endian)
offset 48: magic (uint64_t = 0x4553554D45525F4D)
offset 56: wheel_mod (uint64_t)  — validates resume uses same wheel
```

- `last_sieved` — last number fully processed; resume continues from here + 1.
- `total_count` — cumulative prime count at checkpoint time.
- `entry_count` — number of valid entries in `primes_state.bin`.  On resume,
  the program seeks to `entry_count × 32` before appending, overwriting any
  garbage bytes left by a crashed write.
- `original_target` — target from the initial run.  If a resume target is
  smaller, a warning is printed.
- `magic` — distinguishes a valid checkpoint from arbitrary data
  (`0x4553554D45525F4D` ≙ `"M_RESUME"`).
- `wheel_mod` — wheel modulus used for the sieve.  Resume with a different
  `--wheel` value will error.

## Files

```
fastsieve.c          — source code
README.md            — this file
.gitignore           — ignores binaries, state file, output/
primes_state.bin     — state file (generated, git-ignored)
primes_state.ckpt    — resume checkpoint (generated, git-ignored)
output/              — optional output directory for -o files
```
