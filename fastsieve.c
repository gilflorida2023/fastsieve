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
#include <unistd.h>

/*
 * unsigned __int128 is a GCC/Clang extension providing 128-bit integers.
 * Max safe target: ~3.4e38 (p*p overflow at p > 2^64).
 * Printing requires custom routines (no printf format specifier).
 */
typedef unsigned __int128 u128;

#define FILENAME "primes_state.bin"
#define CHECKPOINT_FILE "primes_state.ckpt"
#define CHECKPOINT_MAGIC ((uint64_t)0x4553554D45525F4D)

static uint64_t WHEEL_MOD = 210;
static uint64_t WHEEL_SIZE = 48;
static uint64_t WHEEL_BYTES = 48;
static uint64_t *wheel_residues = NULL;
static uint16_t *residue_to_bit = NULL;
static uint64_t *wheel_primes = NULL;
static int num_wheel_primes = 0;

/*
 * Bucket sieve: instead of iterating over ALL π(√N) base primes per
 * segment, each prime is placed into one of `nbuckets` circular buckets
 * keyed by which segment its next multiple falls in.  When processing
 * segment S only primes in bucket[S & bucket_mask] are visited — the
 * rest are skipped until their next multiple actually lands in a
 * segment.  For N=10¹² this reduces ~78K outer iterations/segment to
 * ~49, a >1600× reduction.
 *
 * Tomás Oliveira e Silva (2001), also used by Kim Walisch's primesieve.
 */
typedef struct BucketNode {
    u128 prime;
    u128 next_multiple;
    struct BucketNode *next;
} BucketNode;

static BucketNode **buckets = NULL;
static BucketNode *node_pool = NULL;
static uint64_t node_pool_used = 0;
static uint64_t nbuckets = 0;
static uint64_t bucket_mask = 0;

/* Arena allocator: single contiguous region, bump-pointer. No free. */
typedef struct {
    char *base;
    char *ptr;
    char *end;
} Arena;
static Arena aren = {0};

static void *arena_alloc(Arena *a, size_t sz) {
    /* align to 16 bytes */
    sz = (sz + 15) & ~(size_t)15;
    if (a->ptr + sz > a->end) return NULL;
    void *p = a->ptr;
    a->ptr += sz;
    return p;
}
static void arena_reset(void) { if (aren.base) aren.ptr = aren.base; }

static u128 bucket_offset = 0;      /* first number of first segment */
static u128 bucket_bufsize = 0;     /* copy of buffer_size */
static u128 bucket_target = 0;      /* overall target (for drop check) */

static void generate_wheel(uint64_t mod) {
    static const uint64_t valid_wheels[] = {2, 6, 30, 210, 2310};
    int valid = 0;
    for (size_t i = 0; i < sizeof(valid_wheels)/sizeof(valid_wheels[0]); i++) {
        if (valid_wheels[i] == mod) { valid = 1; break; }
    }
    if (!valid) {
        fprintf(stderr, "Invalid wheel modulus: %llu (must be 2, 6, 30, 210, or 2310)\n", mod);
        exit(1);
    }

    WHEEL_MOD = mod;

    uint64_t tmp = mod;
    num_wheel_primes = 0;
    for (uint64_t p = 2; p * p <= tmp; p++) {
        if (tmp % p == 0) {
            num_wheel_primes++;
            while (tmp % p == 0) tmp /= p;
        }
    }
    if (tmp > 1) num_wheel_primes++;

    wheel_primes = malloc(num_wheel_primes * sizeof(uint64_t));
    if (!wheel_primes) { fprintf(stderr, "Out of memory\n"); exit(1); }

    tmp = mod;
    int idx = 0;
    for (uint64_t p = 2; p * p <= tmp; p++) {
        if (tmp % p == 0) {
            wheel_primes[idx++] = p;
            while (tmp % p == 0) tmp /= p;
        }
    }
    if (tmp > 1) wheel_primes[idx++] = tmp;

    WHEEL_SIZE = 0;
    for (uint64_t i = 1; i < mod; i++) {
        uint64_t g = i, t = mod;
        while (t) { uint64_t r = g % t; g = t; t = r; }
        if (g == 1) WHEEL_SIZE++;
    }

    wheel_residues = malloc(WHEEL_SIZE * sizeof(uint64_t));
    residue_to_bit = malloc(WHEEL_MOD * sizeof(uint16_t));
    if (!wheel_residues || !residue_to_bit) { fprintf(stderr, "Out of memory\n"); exit(1); }

    WHEEL_BYTES = WHEEL_SIZE;

    idx = 0;
    for (uint64_t i = 1; i < mod; i++) {
        uint64_t g = i, t = mod;
        while (t) { uint64_t r = g % t; g = t; t = r; }
        if (g == 1) wheel_residues[idx++] = i;
    }

    for (uint64_t i = 0; i < WHEEL_MOD; i++) residue_to_bit[i] = 0xFFFF;
    for (uint64_t i = 0; i < WHEEL_SIZE; i++) residue_to_bit[wheel_residues[i]] = (uint16_t)i;
}

/* --- Bucket sieve helpers ---------------------------------------------- */

/* Return bucket index for a given multiple (low 20 bits of segment num). */
static uint64_t bucket_for(u128 m) {
    u128 seg = (m - bucket_offset) / bucket_bufsize;
    return (uint64_t)(seg & bucket_mask);
}

/* Bump-allocate one node from the pre-allocated pool. */
static BucketNode *node_alloc(void) {
    return &node_pool[node_pool_used++];
}

static void bucket_add(u128 p, u128 nm) {
    BucketNode *n = node_alloc();
    n->prime = p;
    n->next_multiple = nm;
    uint64_t bi = bucket_for(nm);
    n->next = buckets[bi];
    buckets[bi] = n;
}

/* (Re-)initialize the bucket system after bucket_offset/bufsize/target
 * are set.  Allocates the bucket array; primes are added on-the-fly
 * during extend_base_primes and process_segment. */
static bool use_buckets = false;  /* only for large N (threshold in main) */

static void bucket_init(void) {
    node_pool_used = 0;
    memset(buckets, 0, nbuckets * sizeof(BucketNode *));
}

/* --- PrimeEntry --------------------------------------------------------
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

    /* wheel_mod — validates resume uses same wheel */
    uint64_t wheel = WHEEL_MOD;
    fwrite(&wheel, 8, 1, fp);

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

    /* wheel_mod — validate resume uses same wheel */
    uint64_t saved_wheel;
    if (fread(&saved_wheel, 8, 1, fp) != 1) { fclose(fp); return 0; }
    if (saved_wheel != WHEEL_MOD) {
        fprintf(stderr, "Checkpoint wheel mismatch: saved=%llu, current=%llu\n",
                saved_wheel, WHEEL_MOD);
        fclose(fp);
        return 0;
    }

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
    base[base_count].prime = p;
    base[base_count].next = p * 2;
    base_count++;
}

/* Upper bound for π(x) using Rosser's theorem: π(x) < x / (ln x - 4)
 * for x >= 55.  We use a binary-log approximation for ln(x) and add a
 * 30 % safety margin to guarantee the allocation is large enough. */
static uint64_t estimate_pi_max(u128 n) {
    if (n < 55) return 100;
    int bits = 0;
    u128 tmp = n;
    while (tmp) { bits++; tmp >>= 1; }
    double ln_n = (bits - 1) * 0.6931471805599453;
    double est = (double)n / (ln_n - 4.0);
    if (est < 100) return 100;
    uint64_t ret = (uint64_t)(est * 1.3) + 100000;
    if (ret > 2000000000ULL) ret = 2000000000ULL;
    return ret;
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
        for (int i = 0; i < num_wheel_primes; i++) {
            add_base_prime((u128)wheel_primes[i]);
        }
        base_sieved = wheel_primes[num_wheel_primes - 1];
        if (limit <= base_sieved) return;
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
                uint16_t ri = residue_to_bit[r];
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
                    uint16_t ri2 = residue_to_bit[r];
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

                /* Bucket new prime for the main sieve */
                if (use_buckets && buckets && n > (u128)wheel_primes[num_wheel_primes - 1]) {
                    u128 next_m = n * n;
                    if (next_m < bucket_offset) {
                        next_m = bucket_offset;
                        if (next_m % n != 0)
                            next_m += n - next_m % n;
                    }
                    if (next_m <= bucket_target)
                        bucket_add(n, next_m);
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

    if (!use_buckets) {
        /* Linear scan: iterate all base primes for every block */
        for (uint64_t bi = 0; bi < base_count; bi++) {
            u128 p = base[bi].prime;
            u128 m = base[bi].next;
            if (m > seg_end) continue;
            uint64_t step = (uint64_t)(p % WHEEL_MOD);
            u128 blk_step = p / WHEEL_MOD;
            uint64_t r = (uint64_t)(m % WHEEL_MOD);
            u128 block = m / WHEEL_MOD;
            while (m <= seg_end) {
                uint16_t ri = residue_to_bit[r];
                if (ri < WHEEL_SIZE) {
                    uint64_t b = (uint64_t)(block - first_block);
                    if (b < num_blocks)
                        buf[b * WHEEL_BYTES + ri] = 0;
                }
                m += p;
                block += blk_step;
                r += step;
                if (r >= WHEEL_MOD) { r -= WHEEL_MOD; block++; }
            }
            base[bi].next = m;
        }
    } else {
        /* Bucket-driven Phase 1: process primes whose next_multiple
         * falls in this segment, then re-bucket. */
        uint64_t cur_bucket = ((seg_start - bucket_offset) / bucket_bufsize) & bucket_mask;
        BucketNode *node = buckets[cur_bucket];
        buckets[cur_bucket] = NULL;

        while (node) {
            BucketNode *next = node->next;
            u128 p = node->prime;
            u128 m = node->next_multiple;

            if (m <= seg_end) {
                uint64_t step = (uint64_t)(p % WHEEL_MOD);
                u128 blk_step = p / WHEEL_MOD;
                uint64_t r = (uint64_t)(m % WHEEL_MOD);
                u128 block = m / WHEEL_MOD;
                while (m <= seg_end) {
                    uint16_t ri = residue_to_bit[r];
                    if (ri < WHEEL_SIZE) {
                        uint64_t b = (uint64_t)(block - first_block);
                        buf[b * WHEEL_BYTES + ri] = 0;
                    }
                    m += p;
                    block += blk_step;
                    r += step;
                    if (r >= WHEEL_MOD) { r -= WHEEL_MOD; block++; }
                }
                node->next_multiple = m;
            }

            u128 nm = node->next_multiple;
            if (nm <= bucket_target) {
                uint64_t bi = ((nm - bucket_offset) / bucket_bufsize) & bucket_mask;
                node->next = buckets[bi];
                buckets[bi] = node;
            }

            node = next;
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

            u128 next_m;
            if (start_m <= seg_end) {
                u128 m = start_m;
                uint64_t step = (uint64_t)(n % WHEEL_MOD);
                u128 blk_step = n / WHEEL_MOD;
                uint64_t r = (uint64_t)(m % WHEEL_MOD);
                u128 block = m / WHEEL_MOD;
                while (m <= seg_end) {
                    uint16_t ri2 = residue_to_bit[r];
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
                next_m = m;
            } else {
                next_m = start_m;
            }

            if (use_buckets && next_m <= bucket_target)
                bucket_add(n, next_m);
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
    u128 buf = ((u128)262144 / WHEEL_BYTES) * WHEEL_MOD;
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
    printf("Usage: fastsieve [--wheel N] [buffer_size] [target] [-c] [-s] [-r] [-R] [-o file]\n");
    printf("  --wheel N    wheel modulus: 2, 6, 30, 210 (default), 2310\n");
    printf("  buffer_size  segment window size in natural numbers (optional, auto-opt)\n");
    printf("  target       sieve up to this number (required for sieve, report, or resume)\n");
    printf("  -c           count only (no state file, faster)\n");
    printf("  -s           save state file (overrides -c, batched writes)\n");
    printf("  -r           report from existing primes_state.bin (no sieve)\n");
    printf("  -R           resume from checkpoint (requires target > last sieved)\n");
    printf("  -o file      write primes to file (use \"-\" for stdout)\n");
}

int main(int argc, char **argv) {
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--wheel") == 0 && i + 1 < argc) {
            u128 val;
            if (!parse_u128(argv[i + 1], &val)) {
                fprintf(stderr, "Invalid wheel value: %s\n", argv[i + 1]);
                return 1;
            }
            WHEEL_MOD = (uint64_t)val;
        }
    }
    if (WHEEL_MOD != 2 && WHEEL_MOD != 6 && WHEEL_MOD != 30 && WHEEL_MOD != 210 && WHEEL_MOD != 2310) {
        fprintf(stderr, "Invalid wheel: %llu (must be 2, 6, 30, 210, or 2310)\n", WHEEL_MOD);
        return 1;
    }
    generate_wheel(WHEEL_MOD);

    u128 target = 0;
    u128 buffer_size = 0;
    char *output_filename = NULL;
    int count_only = 0;
    int report_mode = 0;
    int save_state = 0;
    int resume_mode = 0;

    int num_count = 0;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--wheel") == 0) {
            i++;
            continue;
        }
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

    /* ----- Arena / memory pre-allocation ---------------------------------- */
    {
        u128 sqrt_target = sqrt_u128(target);
        uint64_t max_base = estimate_pi_max(sqrt_target);
        uint64_t max_nodes = max_base;

        u128 est_segments = (target / buffer_size) + 2;
        nbuckets = 1 << 20;
        while (nbuckets < est_segments) {
            if (nbuckets >= (1ULL << 26)) break;
            nbuckets <<= 1;
        }
        bucket_mask = nbuckets - 1;

        size_t arena_size =
            max_base * sizeof(PrimeEntry) +
            max_nodes * sizeof(BucketNode) +
            nbuckets * sizeof(BucketNode *) +
            FLUSH_BATCH * sizeof(PrimeEntry) +
            256 * 1024 * 1024;

        long pages = sysconf(_SC_PHYS_PAGES);
        long pgsz = sysconf(_SC_PAGESIZE);
        uint64_t avail = (uint64_t)pages * (uint64_t)pgsz / 2;
        if (arena_size > avail) {
            fprintf(stderr, "Target requires ~%zu bytes but only %llu available.\n",
                    arena_size, (unsigned long long)avail);
            return 1;
        }

        char *arena_base = malloc(arena_size);
        if (!arena_base) {
            fprintf(stderr, "Out of memory\n");
            return 1;
        }
        aren.base = arena_base;
        aren.ptr = arena_base;
        aren.end = arena_base + arena_size;

        /* Aligned alloc for cache-line friendly access (64-byte alignment) */
        base = aligned_alloc(64, max_base * sizeof(PrimeEntry));
        if (!base) { fprintf(stderr, "Out of memory\n"); return 1; }
        node_pool = arena_alloc(&aren, max_nodes * sizeof(BucketNode));
        buckets = arena_alloc(&aren, nbuckets * sizeof(BucketNode *));
        write_buf = arena_alloc(&aren, FLUSH_BATCH * sizeof(PrimeEntry));
        memset(buckets, 0, nbuckets * sizeof(BucketNode *));
    }

    FILE *state_fp = NULL;
    if (!count_only || save_state) {
        if (resume_mode) {
            state_fp = fopen(FILENAME, "rb+");
            if (!state_fp) {
                perror("Failed to open " FILENAME " for resume");
                goto fail;
            }
            fseek(state_fp, state_entry_count * 32, SEEK_SET);
        } else {
            state_fp = fopen(FILENAME, "wb+");
            if (!state_fp) {
                perror("Failed to open state file");
                goto fail;
            }
        }
    }

    FILE *output_fp = NULL;
    if (output_filename) {
        if (strcmp(output_filename, "-") == 0) {
            output_fp = stdout;
        } else {
            output_fp = fopen(output_filename, "w");
            if (!output_fp) {
                perror("Failed to open output file");
                if (state_fp) fclose(state_fp);
                goto fail;
            }
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

    /* Wheel primes define the wheel and are NOT in the wheel
     * residues (they'd be excluded as wheel divisors). Count and
     * write them explicitly before the segmented loop begins.
     * Skipped in resume mode — they are already in the state file. */
    u128 first_segment_start = 2;
    u128 last_wheel_prime = 0;
    if (!resume_mode) {
        for (int i = 0; i < num_wheel_primes; i++) {
            if (wheel_primes[i] > target) continue;
            total_count++;
            last_wheel_prime = wheel_primes[i];
            if (state_fp) {
                PrimeEntry e = {wheel_primes[i], (u128)wheel_primes[i] * 2};
                buffer_entry(state_fp, &e);
            }
            if (output_fp) {
                print_u128_f(output_fp, (u128)wheel_primes[i]);
                fputc('\n', output_fp);
            }
        }
        first_segment_start = last_wheel_prime + 1;
        /* Initial checkpoint right after wheel primes */
        if (state_fp) {
            flush_entries(state_fp);
            checkpoint_write(last_wheel_prime, original_target, total_count,
                             state_entry_count);
        }
    } else {
        first_segment_start = current;
    }
    /* Initialize base primes up to sqrt(target) before main sieve */
    extend_base_primes(sqrt_u128(target));
    /* Start main sieve after the wheel primes (not base_sieved) - only for fresh run */
    if (!resume_mode) {
        current = (u128)wheel_primes[num_wheel_primes - 1] + 1;
    }

    /* Initialize bucket sieve (only for targets > 10^11) */
    use_buckets = (target > 100000000000ULL);
    if (use_buckets) {
        u128 first_after_wheel = (u128)wheel_primes[num_wheel_primes - 1] + 1;
        bucket_offset = resume_mode ? res_last + 1 : first_after_wheel;
        bucket_bufsize = buffer_size;
        bucket_target = target;
        bucket_init();
    }

    while (current <= target) {
        u128 seg_end = current + buffer_size - 1;
        if (seg_end > target) seg_end = target;

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
    if (output_fp) fclose(output_fp);

    free(base);
    free(aren.base);
    aren.base = NULL;
    aren.ptr = NULL;
    aren.end = NULL;;
    return 0;

fail:
    free(base);
    free(aren.base);
    arena_reset();
    return 1;
}
