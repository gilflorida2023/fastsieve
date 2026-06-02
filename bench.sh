#!/usr/bin/env bash
set -euo pipefail

RESULTS="bench_results.txt"
SRC="fastsieve.c"
BIN="./fastsieve"
CKPT="primes_state.ckpt"

rm -f "$RESULTS"
echo "fastsieve benchmark report" | tee "$RESULTS"
echo "=========================" | tee -a "$RESULTS"
echo "host: $(uname -a)" | tee -a "$RESULTS"
echo "gcc:  $(gcc --version | head -1)" | tee -a "$RESULTS"
echo "buffer: auto_opt_buffer() => ~$(./fastsieve 100000000000 -c 2>&1 | grep Buffer | awk '{print $2}')" | tee -a "$RESULTS"
echo "" | tee -a "$RESULTS"

# ── 1. Build ──────────────────────────────────────────────────────────────
echo "── Build (O3, fastest per benchmarks) ──" | tee -a "$RESULTS"
gcc -O3 -pipe -Wall -Wextra -o $BIN $SRC 2>&1
echo "  OK" | tee -a "$RESULTS"
echo "" | tee -a "$RESULTS"

# ── 2. Accuracy ──────────────────────────────────────────────────────────
echo "── Accuracy ──" | tee -a "$RESULTS"
EXPECT=("1000:168" "1000000:78498" "10000000:664579" "100000000:5761455" "1000000000:50847534")
ok=1
for pair in "${EXPECT[@]}"; do
    t="${pair%%:*}"
    expected="${pair##*:}"
    got=$($BIN "$t" -c 2>&1 | sed -n 's/Found \([0-9]*\) primes.*/\1/p')
    if [[ "$got" != "$expected" ]]; then
        echo "  \u03c0($t): FAIL (got $got, expected $expected)" | tee -a "$RESULTS"
        ok=0
    else
        echo "  \u03c0($t) = $got PASS" | tee -a "$RESULTS"
    fi
done
[[ "$ok" == "1" ]] && echo "  all pass" | tee -a "$RESULTS"
echo "" | tee -a "$RESULTS"

# ── 3. Fresh speed runs (no state needed) ──────────────────────────────────
echo "── Speed (fresh -c, 5 runs each) ──" | tee -a "$RESULTS"

bench_one() {
    local target="$1" runs="$2"
    local best=99999 sum=0 val
    for i in $(seq 1 "$runs"); do
        val=$(/usr/bin/time -f "%e" $BIN "$target" -c 2>&1 >/dev/null)
        sum=$(echo "$sum + $val" | bc)
        if (( $(echo "$val < $best" | bc) )); then best=$val; fi
    done
    avg=$(echo "scale=4; $sum / $runs" | bc)
    echo "$target min=$best avg=$avg" >> "$RESULTS"
    printf "  %-14s min=%.3f  avg=%.3f\n" "$target" "$best" "$avg"
}

for t in 1000000 10000000 100000000 1000000000 10000000000; do
    bench_one "$t" 5
done
echo "" | tee -a "$RESULTS"

# ── 4. Long runs ────────────────────────────────────────────────────────
echo "── Long runs (no state file: -c only) ──" | tee -a "$RESULTS"

do_run() {
    local target="$1" label="$2"
    echo "  $label ..." | tee -a "$RESULTS"
    set +e
    output=$($BIN -c "$target" 2>&1)
    rc=$?
    set -e
    if [[ $rc -eq 0 ]]; then
        dur=$(echo "$output" | sed -n 's/.*Found.*(\(.*\))/\1/p')
        echo "  $label done => $dur" | tee -a "$RESULTS"
        echo "$label target=$target dur=$dur" >> "$RESULTS"
    else
        echo "  $label FAILED (exit $rc)" | tee -a "$RESULTS"
        echo "$label target=$target result=fail,exit=$rc" >> "$RESULTS"
    fi
}

do_run "100000000000" "10^11"
do_run "1000000000000" "10^12"
echo "" | tee -a "$RESULTS"

echo "Done. Results in $RESULTS" | tee -a "$RESULTS"
