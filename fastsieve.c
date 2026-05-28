/*
 * fastsieve — Segmented Sieve of Eratosthenes with mod-210 Wheel
 *
 * Finds/counts primes in ranges up to ~10^38 using:
 *   - mod-210 wheel factorization: pre-filters numbers divisible by
 *     2,3,5,7, leaving only 48 of 210 candidates (77% reduction)
 *   - Segmented processing: fixed-size windows keep the sieve buffer
 *     L1-cache-friendly (~32 KB), requiring only primes up to sqrt(N)
 *     in memory at once
 *   - unsigned __int128: extends safe range beyond what uint64_t
 *     allows (p*p overflow at p > 2^32)
 *   - Dynamic base primes: bootstraps from {2,3,5,7} via squaring,
 *     growing the base-prime list as sqrt(segment_end) increases
 *
 * Compile:  gcc -O3 -march=native -o fastsieve fastsieve.c
 * Example:  ./fastsieve -c 1000000000000
 *           ./fastsieve -r
 *           ./fastsieve -c -o primes.txt 1000000000
 *
 * Reference: Pritchard (1981) "A Sublinear Additive Sieve for
 * Finding Prime Numbers" — wheel factorization origin.
 * Modern practice: Kim Walisch's primesieve, Tomás Oliveira e Silva.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>

/*
 * unsigned __int128 is a GCC/Clang extension providing 128-bit integers.
 * Max safe target: ~3.4e38 (p*p overflow at p > 2^64).
 * Printing requires custom routines (no printf format specifier).
 */
typedef unsigned __int128 u128;

#define FILENAME "primes_state.bin"
#define WHEEL_MOD 210
#define WHEEL_SIZE 48
#define WHEEL_BYTES 48

/*
 * mod-210 wheel residues: numbers < 210 that are coprime to 210.
 * 210 = 2 * 3 * 5 * 7, so these 48 residues are not divisible by
 * 2, 3, 5, or 7 — the only possible prime values modulo 210.
 * Every prime >= 11 falls into exactly one of these residues.
 * This eliminates 77% of candidates before sieving begins.
 */
static const uint64_t wheel_residues[WHEEL_SIZE] = {
    1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
    53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 121, 127, 131,
    137, 139, 143, 149, 151, 157, 163, 167, 169, 173,
    179, 181, 187, 191, 193, 197, 199, 209
};

static uint8_t residue_to_bit[WHEEL_MOD];

static void build_lut(void) {
    for (int i = 0; i < WHEEL_MOD; i++)
        residue_to_bit[i] = 0xFF;
    for (int i = 0; i < WHEEL_SIZE; i++)
        residue_to_bit[wheel_residues[i]] = i;
}

/*
 * PrimeEntry: one entry in the binary state file (primes_state.bin).
 *   prime: a discovered prime value
 *   next:  first multiple of prime that falls in or after the
 *          current segment (used for resumability)
 *
 * Layout on disk: 4 consecutive uint64_t values in little-endian:
 *   [prime_lo, prime_hi, next_lo, next_hi]  (32 bytes total)
 */
typedef struct {
    u128 prime;
    u128 next;
} PrimeEntry;

static PrimeEntry *base = NULL;
static uint64_t base_count = 0;
static uint64_t base_cap = 0;
static u128 base_sieved = 0;

static uint64_t total_count = 0;

/* Batched state file writes: accumulate PrimeEntry structs in memory
 * and flush to disk once per segment (every ~65K primes), reducing
 * syscall overhead vs one fwrite per prime. */
#define FLUSH_BATCH 65536
static PrimeEntry *write_buf = NULL;
static uint64_t write_buf_count = 0;

/* Track how many entries have been flushed to the state file.
 * Used by the resume checkpoint to know where valid entries end
 * (entries beyond this count may be garbage from a crashed write). */
static uint64_t state_entry_count = 0;

/* -- Resume checkpoint support --------------------------------------------
 * A separate small file (primes_state.ckpt) stores the last committed
 * resume position.  It is atomically updated via write-to-temp + rename,
 * so a crash never corrupts the checkpoint (the previous copy survives).
 *
 * Layout on disk (48 bytes total):
 *   [last_sieved_lo (8)] [last_sieved_hi (8)]   -- u128
 *   [total_count (8)]                            -- uint64_t
 *   [entry_count (8)]                            -- uint64_t
 *   [original_target_lo (8)] [original_target_hi (8)]  -- u128
 *   [magic (8)]                                  -- uint64_t
 */
#define CHECKPOINT_FILE "primes_state.ckpt"
#define CHECKPOINT_MAGIC ((uint64_t)0x4553554D45525F4D) /* "M_RESUME" */

static void print_u128(u128 n) {
    if (n == 0) { putchar('0'); return; }
    char buf[48];
    int i = 47;
    buf[i] = '\0';
    while (n > 0) {
        buf[--i] = '0' + (char)(n % 10);
        n /= 10;
    }
    fputs(buf + i, stdout);
}

static void print_u128_f(FILE *fp, u128 n) {
    if (n == 0) { fputc('0', fp); return; }
    char buf[48];
    int i = 47;
    buf[i] = '\0';
    while (n > 0) {
        buf[--i] = '0' + (char)(n % 10);
        n /= 10;
    }
    fputs(buf + i, fp);
}

static u128 sqrt_u128(u128 n) {
    if (n < 2) return n;
    u128 lo = 1, hi = (u128)1 << 64;
    if (hi > n) hi = n;
    while (lo < hi) {
        u128 mid = lo + (hi - lo) / 2;
        if (mid <= n / mid && (mid + 1) > n / (mid + 1))
            return mid;
        if (mid > n / mid)
            hi = mid;
        else
            lo = mid + 1;
    }
    return lo;
}

static int parse_u128(const char *s, u128 *out) {
    *out = 0;
    if (!s || !*s) return 0;
    while (*s) {
        if (*s < '0' || *s > '9') return 0;
        *out = *out * 10 + (u128)(*s - '0');
        s++;
    }
    return 1;
}

/* Flush buffered entries to disk (portable serialization: each u128
 * as two uint64_t halves in little-endian). Must be called before
 * fclose and before any fseek/read on the state file. */
static void flush_entries(FILE *state_fp) {
    if (!state_fp || write_buf_count == 0) return;
    for (uint64_t i = 0; i < write_buf_count; i++) {
        PrimeEntry *e = &write_buf[i];
        uint64_t lo, hi;
        lo = (uint64_t)e->prime;
        hi = (uint64_t)(e->prime >> 64);
        fwrite(&lo, 8, 1, state_fp);
        fwrite(&hi, 8, 1, state_fp);
        lo = (uint64_t)e->next;
        hi = (uint64_t)(e->next >> 64);
        fwrite(&lo, 8, 1, state_fp);
        fwrite(&hi, 8, 1, state_fp);
    }
    state_entry_count += write_buf_count;
    write_buf_count = 0;
}

/* Buffer a PrimeEntry for batched writing. Flushes to disk when the
 * buffer fills. state_fp may be NULL (no-op if no state file). */
static void buffer_entry(FILE *state_fp, PrimeEntry *e) {
    if (!state_fp) return;
    write_buf[write_buf_count++] = *e;
    if (write_buf_count >= FLUSH_BATCH)
        flush_entries(state_fp);
}

static void checkpoint_write(u128 last_sieved, u128 original_target,
                             uint64_t total_count, uint64_t entry_count) {
    FILE *fp = fopen(CHECKPOINT_FILE ".tmp", "wb");
    if (!fp) return;

    /* serialize last_sieved (u128 as two uint64_t halves) */
    uint64_t lo = (uint64_t)last_sieved;
    uint64_t hi = (uint64_t)(last_sieved >> 64);
    fwrite(&lo, 8, 1, fp);
    fwrite(&hi, 8, 1, fp);

    /* total_count, entry_count */
    fwrite(&total_count, 8, 1, fp);
    fwrite(&entry_count, 8, 1, fp);

    /* original_target (u128) */
    lo = (uint64_t)original_target;
    hi = (uint64_t)(original_target >> 64);
    fwrite(&lo, 8, 1, fp);
    fwrite(&hi, 8, 1, fp);

    /* magic — marks a complete, valid footer */
    uint64_t magic = CHECKPOINT_MAGIC;
    fwrite(&magic, 8, 1, fp);

    fclose(fp);

    /* Atomic replacement: the rename(2) is POSIX-atomic on the same
     * filesystem.  If interrupted here, the original .ckpt still exists
     * and remains valid. */
    rename(CHECKPOINT_FILE ".tmp", CHECKPOINT_FILE);
}

static int checkpoint_read(u128 *last_sieved, u128 *original_target,
                           uint64_t *total_count, uint64_t *entry_count) {
    FILE *fp = fopen(CHECKPOINT_FILE, "rb");
    if (!fp) return 0;

    uint64_t lo, hi, magic;
    if (fread(&lo, 8, 1, fp) != 1) { fclose(fp); return 0; }
    if (fread(&hi, 8, 1, fp) != 1) { fclose(fp); return 0; }
    *last_sieved = ((u128)hi << 64) | lo;

    if (fread(total_count, 8, 1, fp) != 1) { fclose(fp); return 0; }
    if (fread(entry_count, 8, 1, fp) != 1) { fclose(fp); return 0; }

    if (fread(&lo, 8, 1, fp) != 1) { fclose(fp); return 0; }
    if (fread(&hi, 8, 1, fp) != 1) { fclose(fp); return 0; }
    *original_target = ((u128)hi << 64) | lo;

    if (fread(&magic, 8, 1, fp) != 1) { fclose(fp); return 0; }
    if (magic != CHECKPOINT_MAGIC)    { fclose(fp); return 0; }

    fclose(fp);
    return 1;
}

static int read_entry(FILE *fp, PrimeEntry *e) {
    uint64_t lo, hi;
    if (fread(&lo, 8, 1, fp) != 1) return 0;
    if (fread(&hi, 8, 1, fp) != 1) return 0;
    e->prime = ((u128)hi << 64) | lo;
    if (fread(&lo, 8, 1, fp) != 1) return 0;
    if (fread(&hi, 8, 1, fp) != 1) return 0;
    e->next = ((u128)hi << 64) | lo;
    return 1;
}

static void add_base_prime(u128 p) {
    if (base_count >= base_cap) {
        base_cap = base_cap ? base_cap * 2 : 65536;
        base = realloc(base, base_cap * sizeof(PrimeEntry));
    }
    base[base_count].prime = p;
    base[base_count].next = 0;
    base_count++;
}

/*
 * Grow the base-prime list so it covers all primes up to `limit`.
 * Uses bootstrapping by squaring:
 *   Level 0: hardcode {2,3,5,7}
 *   Level 1: sieve [8, 49]     using level 0 (7^2 = 49)
 *   Level 2: sieve [50, 2401]  using level 1 (49^2 = 2401)
 *   ...
 * This converges rapidly (5 levels exceed 128-bit range).
 * Within each level, the same wheel-segmented approach is used:
 *   - Clear a buffer of wheel candidates
 *   - Mark composites using existing base primes
 *   - Scan survivors: each is a new base prime
 *   - Mark multiples of each new prime in remaining buffer
 */
static void extend_base_primes(u128 limit) {
    if (limit <= base_sieved) return;
    if (limit > (u128)1 << 62) limit = (u128)1 << 62;

    if (base_count == 0) {
        add_base_prime(2);
        add_base_prime(3);
        add_base_prime(5);
        add_base_prime(7);
        base_sieved = 7;
        if (limit <= 7) return;
    }

    while (base_sieved < limit) {
        u128 next_target;
        if (base_sieved >= ((u128)1 << 62))
            next_target = limit;
        else
            next_target = base_sieved * base_sieved;
        if (next_target > limit || next_target <= base_sieved)
            next_target = limit;
        if (next_target <= base_sieved) break;

        u128 start = base_sieved + 1;
        u128 end = next_target;
        if (start > end) break;

        uint64_t first_block = (uint64_t)(start / WHEEL_MOD);
        uint64_t last_block = (uint64_t)(end / WHEEL_MOD);
        uint64_t num_blocks = last_block - first_block + 1;

        uint8_t *buf = malloc(num_blocks * WHEEL_BYTES);
        if (!buf) { fprintf(stderr, "Out of memory\n"); exit(1); }
        memset(buf, 0xFF, num_blocks * WHEEL_BYTES);

        for (uint64_t bi = 0; bi < base_count; bi++) {
            u128 p = base[bi].prime;
            if (p < 11) continue;
            u128 fm = start;
            if (p * p > fm) fm = p * p;
            if (fm > end) continue;
            if (fm % p != 0)
                fm += p - fm % p;
            u128 m = fm;
            uint64_t step = (uint64_t)(p % WHEEL_MOD);
            u128 blk_step = p / WHEEL_MOD;
            uint64_t r = (uint64_t)(m % WHEEL_MOD);
            u128 block = m / WHEEL_MOD;
            while (m <= end) {
                uint8_t ri = residue_to_bit[r];
                if (ri < WHEEL_SIZE) {
                    uint64_t buf_i = (uint64_t)(block - first_block);
                    buf[buf_i * WHEEL_BYTES + ri] = 0;
                }
                m += p;
                block += blk_step;
                r += step;
                if (r >= WHEEL_MOD) { r -= WHEEL_MOD; block++; }
            }
        }

        for (uint64_t bi = 0; bi < num_blocks; bi++) {
            for (int ri = 0; ri < WHEEL_SIZE; ri++) {
                if (!buf[bi * WHEEL_BYTES + ri]) continue;
                u128 n = (first_block + bi) * (u128)WHEEL_MOD + wheel_residues[ri];
                if (n < start || n > end) continue;

                add_base_prime(n);

                u128 start_m = n * n;
                if (start_m < start) start_m = start;
                u128 m = start_m;
                uint64_t step = (uint64_t)(n % WHEEL_MOD);
                u128 blk_step = n / WHEEL_MOD;
                uint64_t r = (uint64_t)(m % WHEEL_MOD);
                u128 block = m / WHEEL_MOD;
                while (m <= end) {
                    uint8_t ri2 = residue_to_bit[r];
                    if (ri2 < WHEEL_SIZE) {
                        uint64_t b = (uint64_t)(block - first_block);
                        if (b < num_blocks)
                            buf[b * WHEEL_BYTES + ri2] = 0;
                    }
                    m += n;
                    block += blk_step;
                    r += step;
                    if (r >= WHEEL_MOD) { r -= WHEEL_MOD; block++; }
                }
            }
        }

        free(buf);
        base_sieved = end;
    }
}

/*
 * Process one segment [seg_start, seg_end] of the number line.
 *
 * Phase 1 — Mark composites:
 *   For each base prime p, walk its wheel-aligned multiples in the
 *   segment and clear their residue slots in the buffer. The starting
 *   multiple is max(p*p, next_multiple_ge_seg_start) to avoid
 *   redundant work (smaller multiples were handled by smaller primes
 *   or in previous segments).
 *
 * Phase 2 — Discover primes:
 *   Scan the buffer for surviving (non-cleared) residue slots.
 *   Each survivor is a new prime. For each:
 *     - Increment the total count
 *     - Write to the state file (if not -c) and/or output file (if -o)
 *     - Mark its multiples in the remaining buffer (starting from p*p)
 *       so later survivors are not fooled by composites.
 *
 * The buffer uses 1 byte per wheel residue position, 48 bytes per
 * 210-number block. The block layout allows O(1) indexing from any
 * number via: block = number / 210, residue = number % 210.
 */
static void process_segment(FILE *state_fp, u128 seg_start, u128 seg_end,
                            FILE *output_fp) {
    if (seg_start > seg_end) return;

    uint64_t first_block = (uint64_t)(seg_start / WHEEL_MOD);
    uint64_t last_block = (uint64_t)(seg_end / WHEEL_MOD);
    uint64_t num_blocks = last_block - first_block + 1;

    uint8_t *buf = malloc(num_blocks * WHEEL_BYTES);
    if (!buf) { fprintf(stderr, "Out of memory\n"); exit(1); }
    memset(buf, 0xFF, num_blocks * WHEEL_BYTES);

    for (uint64_t i = 0; i < base_count; i++) {
        u128 p = base[i].prime;
        u128 fm = seg_start;
        if (p * p > fm) fm = p * p;
        if (fm > seg_end) continue;
        if (fm % p != 0)
            fm += p - fm % p;

        u128 m = fm;
        uint64_t step = (uint64_t)(p % WHEEL_MOD);
        u128 blk_step = p / WHEEL_MOD;
        uint64_t r = (uint64_t)(m % WHEEL_MOD);
        u128 block = m / WHEEL_MOD;
        while (m <= seg_end) {
            uint8_t ri = residue_to_bit[r];
            if (ri < WHEEL_SIZE) {
                uint64_t b = (uint64_t)(block - first_block);
                buf[b * WHEEL_BYTES + ri] = 0;
            }
            m += p;
            block += blk_step;
            r += step;
            if (r >= WHEEL_MOD) { r -= WHEEL_MOD; block++; }
        }
    }

    for (uint64_t bi = 0; bi < num_blocks; bi++) {
        for (int ri = 0; ri < WHEEL_SIZE; ri++) {
            if (!buf[bi * WHEEL_BYTES + ri]) continue;
            u128 n = (first_block + bi) * (u128)WHEEL_MOD + wheel_residues[ri];
            if (n < seg_start || n > seg_end) continue;

            total_count++;

            PrimeEntry e;
            e.prime = n;
            e.next = n * 2;
            while (e.next <= seg_end)
                e.next += n;

            buffer_entry(state_fp, &e);

            if (output_fp) {
                print_u128_f(output_fp, n);
                fputc('\n', output_fp);
            }

            u128 start_m = n * n;
            if (start_m < seg_start) start_m = seg_start;
            if (start_m > seg_end) continue;
            u128 m = start_m;
            uint64_t step = (uint64_t)(n % WHEEL_MOD);
            u128 blk_step = n / WHEEL_MOD;
            uint64_t r = (uint64_t)(m % WHEEL_MOD);
            u128 block = m / WHEEL_MOD;
            while (m <= seg_end) {
                uint8_t ri2 = residue_to_bit[r];
                if (ri2 < WHEEL_SIZE) {
                    uint64_t b = (uint64_t)(block - first_block);
                    if (b < num_blocks)
                        buf[b * WHEEL_BYTES + ri2] = 0;
                }
                m += n;
                block += blk_step;
                r += step;
                if (r >= WHEEL_MOD) { r -= WHEEL_MOD; block++; }
            }
        }
    }

    free(buf);
}

/*
 * Auto-optimize buffer size: target a 32 KiB sieve buffer.
 *
 * The sieve buffer uses 48 bytes per 210-number block (one byte per
 * wheel residue). To fill roughly 32 KiB (32768 bytes) of L1 cache:
 *   blocks  = 32768 / 48  ≈ 682
 *   numbers = 682 × 210  = 143,220
 *
 * This is clamped to [210, target] so single-segment runs work
 * correctly for small targets.
 */
static u128 auto_opt_buffer(u128 target) {
    u128 buf = ((u128)32768 / WHEEL_BYTES) * WHEEL_MOD;
    if (buf < WHEEL_MOD) buf = WHEEL_MOD;
    if (buf > target) buf = target;
    return buf;
}

static void format_duration(double sec, char *buf, size_t sz) {
    if (sec < 60.0) {
        snprintf(buf, sz, "%.4f sec", sec);
    } else if (sec < 3600.0) {
        int m = (int)(sec / 60);
        double s = sec - m * 60;
        snprintf(buf, sz, "%dm %.1fs", m, s);
    } else if (sec < 86400.0) {
        int h = (int)(sec / 3600);
        int m = (int)((sec - h * 3600) / 60);
        int s = (int)(sec - h * 3600 - m * 60);
        snprintf(buf, sz, "%dh %dm %ds", h, m, s);
    } else {
        int d = (int)(sec / 86400);
        int h = (int)((sec - d * 86400) / 3600);
        int m = (int)((sec - d * 86400 - h * 3600) / 60);
        snprintf(buf, sz, "%dd %dh %dm", d, h, m);
    }
}

static void print_help(void) {
    printf("Usage: fastsieve [buffer_size] [target] [-c] [-s] [-r] [-R] [-o file]\n");
    printf("  buffer_size  segment window size in natural numbers (optional, auto-opt)\n");
    printf("  target       sieve up to this number (required for sieve, report, or resume)\n");
    printf("  -c           count only (no state file, faster)\n");
    printf("  -s           save state file (overrides -c, batched writes)\n");
    printf("  -r           report from existing primes_state.bin (no sieve)\n");
    printf("  -R           resume from checkpoint (requires target > last sieved)\n");
    printf("  -o file      write primes to file (use \"-\" for stdout)\n");
}

int main(int argc, char **argv) {
    build_lut();

    u128 target = 0;
    u128 buffer_size = 0;
    char *output_filename = NULL;
    int count_only = 0;
    int report_mode = 0;
    int save_state = 0;
    int resume_mode = 0;

    int num_count = 0;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-c") == 0) {
            count_only = 1;
        } else if (strcmp(argv[i], "-s") == 0) {
            save_state = 1;
        } else if (strcmp(argv[i], "-r") == 0) {
            report_mode = 1;
        } else if (strcmp(argv[i], "-R") == 0) {
            resume_mode = 1;
        } else if (strcmp(argv[i], "-o") == 0) {
            if (i + 1 < argc) output_filename = argv[++i];
            else { fprintf(stderr, "-o requires a filename\n"); return 1; }
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            print_help();
            return 1;
        } else {
            u128 val;
            if (!parse_u128(argv[i], &val)) {
                fprintf(stderr, "Invalid number: %s\n", argv[i]);
                return 1;
            }
            if (num_count == 0) buffer_size = val;
            else if (num_count == 1) target = val;
            else { fprintf(stderr, "Unexpected argument: %s\n", argv[i]); return 1; }
            num_count++;
        }
    }
    if (num_count == 1) { target = buffer_size; buffer_size = 0; }

    if (report_mode && (count_only || save_state || resume_mode)) {
        fprintf(stderr, "-r cannot be combined with -c, -s, or -R\n");
        return 1;
    }

    if (report_mode) {
        FILE *rfp = fopen(FILENAME, "rb");
        if (!rfp) { perror("Failed to open " FILENAME); return 1; }
        fseek(rfp, 0, SEEK_END);
        long fsize = ftell(rfp);
        rewind(rfp);

        PrimeEntry e;
        u128 first = 0, last = 0;
        uint64_t nentries = 0;
        while (read_entry(rfp, &e)) {
            if (nentries == 0) first = e.prime;
            last = e.prime;
            nentries++;
        }

        printf("State file: %s (%ld bytes)\n", FILENAME, fsize);
        printf("Entries:    %llu\n", (unsigned long long)nentries);
        printf("First:      ");
        print_u128(first);
        printf("\nLast:       ");
        print_u128(last);
        printf("\n");

        if (output_filename) {
            FILE *ofp;
            if (strcmp(output_filename, "-") == 0) {
                ofp = stdout;
            } else {
                ofp = fopen(output_filename, "w");
                if (!ofp) { perror(output_filename); return 1; }
            }
            rewind(rfp);
            while (read_entry(rfp, &e)) {
                print_u128_f(ofp, e.prime);
                fputc('\n', ofp);
            }
            if (ofp != stdout) fclose(ofp);
        }
        fclose(rfp);
        return 0;
    }

    u128 original_target = 0;
    u128 res_last = 0;  /* set in resume-mode block, read later */

    /* ----- Resume mode -------------------------------------------------- */
    if (resume_mode) {
        u128 res_orig;
        uint64_t res_total, res_entries;

        if (!checkpoint_read(&res_last, &res_orig, &res_total, &res_entries)) {
            fprintf(stderr, "No valid checkpoint found in %s.\n"
                            "Run without -R for a fresh sieve, or use -r"
                            " to report.\n", CHECKPOINT_FILE);
            return 1;
        }

        if (target <= res_last) {
            printf("Target already fully sieved (last sieved: ");
            print_u128(res_last);
            printf(").\n");
            print_u128(res_total);
            printf(" primes found.\n");
            return 0;
        }

        if (target < res_orig) {
            printf("Note: new target is smaller than the original target (");
            print_u128(res_orig);
            printf("). Continuing.\n");
        }

        total_count = res_total;
        state_entry_count = res_entries;
        original_target = res_orig;
        base_count = 0;
        base_sieved = 0;

        if (buffer_size == 0)
            buffer_size = auto_opt_buffer(target);
    } else {
        original_target = target;
    }

    if (target == 0) {
        print_help();
        return 1;
    }

    if (buffer_size == 0)
        buffer_size = auto_opt_buffer(target);

    printf("Target: ");
    print_u128(target);
    printf("\nBuffer: ");
    print_u128(buffer_size);
    printf("\n");

    FILE *state_fp = NULL;
    if (!count_only || save_state) {
        if (resume_mode) {
            state_fp = fopen(FILENAME, "rb+");
            if (!state_fp) {
                perror("Failed to open " FILENAME " for resume");
                return 1;
            }
            fseek(state_fp, state_entry_count * 32, SEEK_SET);
        } else {
            state_fp = fopen(FILENAME, "wb+");
            if (!state_fp) {
                perror("Failed to open state file");
                return 1;
            }
        }
        write_buf = malloc(FLUSH_BATCH * sizeof(PrimeEntry));
        if (!write_buf) {
            fprintf(stderr, "Out of memory\n");
            if (state_fp) fclose(state_fp);
            return 1;
        }
    }

    FILE *output_fp = NULL;
    if (output_filename) {
        output_fp = fopen(output_filename, "w");
        if (!output_fp) {
            perror("Failed to open output file");
            if (state_fp) fclose(state_fp);
            return 1;
        }
    }

    clock_t start_time = clock();

    u128 current = 2;
    u128 sqrt_cur = 2;
    uint64_t segs = 0;

    if (resume_mode) {
        current = res_last + 1;
        printf("Resuming from ");
        print_u128(current);
        printf(" to ");
        print_u128(target);
        printf("\n");
    }

    /* Primes 2,3,5,7 define the wheel and are NOT in the wheel
     * residues (they'd be excluded as wheel divisors). Count and
     * write them explicitly before the segmented loop begins.
     * Skipped in resume mode — they are already in the state file. */
    if (!resume_mode) {
        static const u128 small_primes[] = {2, 3, 5, 7};
        for (int i = 0; i < 4; i++) {
            if (small_primes[i] > target) continue;
            total_count++;
            if (state_fp) {
                PrimeEntry e = {small_primes[i], small_primes[i] * 2};
                buffer_entry(state_fp, &e);
            }
            if (output_fp) {
                print_u128_f(output_fp, small_primes[i]);
                fputc('\n', output_fp);
            }
        }
        /* Initial checkpoint right after small primes */
        if (state_fp) {
            flush_entries(state_fp);
            checkpoint_write(7, original_target, total_count,
                             state_entry_count);
        }
    }

    while (current <= target) {
        u128 seg_end = current + buffer_size - 1;
        if (seg_end > target) seg_end = target;

        sqrt_cur = sqrt_u128(seg_end);
        if (sqrt_cur > base_sieved)
            extend_base_primes(sqrt_cur);

        process_segment(state_fp, current, seg_end, output_fp);

        /* Periodic checkpoint (every 10 segments) so a crash loses
         * at most 10 segments of progress. */
        if (state_fp && (segs % 10 == 0)) {
            flush_entries(state_fp);
            checkpoint_write(seg_end, original_target, total_count,
                             state_entry_count);
        }

        current = seg_end + 1;
        segs++;
    }

    clock_t end_time = clock();
    double elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    printf("Found ");
    print_u128(total_count);
    char dur[64];
    format_duration(elapsed, dur, sizeof(dur));
    printf(" primes in %llu segments (%s)\n",
           (unsigned long long)segs, dur);

    if (state_fp) {
        flush_entries(state_fp);
        checkpoint_write(target, original_target, total_count,
                         state_entry_count);
        fclose(state_fp);
    }
    free(write_buf);
    if (output_fp) fclose(output_fp);

    free(base);
    return 0;
}
